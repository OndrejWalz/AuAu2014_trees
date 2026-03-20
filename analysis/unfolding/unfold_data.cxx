// unfold_data.C
#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldSvd.h"
#include "TSVDUnfold.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TMatrixD.h"

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

using std::string;
using std::vector;
using std::cout;
using std::endl;

// ------------------------- config -----------------------------

// Same iterations as closure
static const int kBayesIterDefault = 4;

// SVD regularization values (same as embedding code)
static const int kSVDRegValues[] = {4};
static const int kNSVDReg        = sizeof(kSVDRegValues)/sizeof(kSVDRegValues[0]);

// measured & truth binning (must be IDENTICAL to what you used for response)
static const int nbins_meas = 24;
static const double bin_meas_edges[nbins_meas+1] = {
  -100,-80,-60,-40,-20,-10,-5,-2.5,0,2.5,5,7.5,10,12.5,15,17.5,
  20,22.5,25,27.5,30,35,40,50,60
};

static const int nbins_truth = 10;
static const double bin_truth_edges[nbins_truth+1] = {
  0,5,10,15,20,25,30,35,40,50,60
};

// same centrality & R labels as before
static const vector<string> kCentralities =
  {"CENT_0_10", "MID_20_40", "PERI_60_80"};
static const vector<string> kRadii =
  {"R0.2", "R0.3", "R0.4"};

// pt_lead cuts (same as in closure)
static const double kPtLeadCuts[] = {0.0, 5.0, 7.0, 9.0};
static const int    kNPtLeadCuts  = sizeof(kPtLeadCuts)/sizeof(kPtLeadCuts[0]);

// ------------------ jet-quality cuts -------------------------

static const double CUT_AREA_02 = 0.07;  // R = 0.2
static const double CUT_AREA_03 = 0.20;  // R = 0.3
static const double CUT_AREA_04 = 0.40;  // R = 0.4
static const double CUT_NEUTRAL_FRACTION = 0.95;

// --------------------------------------------------------------

static void EnsureDir(const string& path){
  if (gSystem->AccessPathName(path.c_str()))
    gSystem->mkdir(path.c_str(), /*recursive=*/true);
}

static bool SameBinning1D(const TH1* a, const TH1* b, double tol = 1e-6);

static TH1D* GetEffHistChecked(TFile* fEff,
                               const std::string& tagdirName,
                               const char* objName,
                               const TH1* ref,
                               const char* cloneNameForOut,
                               bool verbose = true);

static void ApplyDivideCorr(TH1D* spec,
                            const TH1D* eff,
                            bool verbose = true);

static void ApplyMultiplyCorr(TH1D* spec,
                              const TH1D* fac,
                              bool verbose = true);                            

/**
 * dataFile     : real data jets (with JetTree)
 * responseFile : single ROOT file containing per-tag directories:
 *                R*_CENT_*_ptlead* / { response, hPrior, ... }
 * outDir       : where unfolded data ROOT file will be written
 * nIter        : Bayes iterations (default 4)
 */
void unfold_data(const char* dataFile,
                 const char* responseFile,
                 const char* effFile,
                 const char* outDir,
                 const char* method = "BAYES",
                 int nIter = kBayesIterDefault)
{

  std::string m(method);
  std::cout << ">>> Unfolding method = " << m << std::endl;


  gStyle->SetOptStat(0);
  EnsureDir(outDir);
  TH1::SetDefaultSumw2(kTRUE);

  TFile* fData = TFile::Open(dataFile, "READ");
  if (!fData || fData->IsZombie()) {
    cout << "[error] Cannot open data file: " << dataFile << endl;
    return;
  }

  TFile* fRespAll = TFile::Open(responseFile, "READ");
  if (!fRespAll || fRespAll->IsZombie()) {
    cout << "[error] Cannot open response file: " << responseFile << endl;
    fData->Close(); delete fData;
    return;
  }

  TFile* fEff = TFile::Open(effFile, "READ");
  if (!fEff || fEff->IsZombie()) {
    cout << "[error] Cannot open efficiencies file: " << effFile << endl;
    fRespAll->Close(); delete fRespAll;
    fData->Close();    delete fData;
    return;
  }

  const string outRootPath = string(outDir) + "/unfolded_data.root";
  TFile* fOutAll = TFile::Open(outRootPath.c_str(), "RECREATE");
  if (!fOutAll || fOutAll->IsZombie()) {
    cout << "[error] Cannot create output file: " << outRootPath << endl;
    fEff->Close();     delete fEff;
    fRespAll->Close(); delete fRespAll;
    fData->Close();    delete fData;
    return;
  }

  cout << "\n=== Unfolding real data ===\n";
  cout << "  Data file     : " << dataFile << "\n";
  cout << "  Response file : " << responseFile << "\n";
  cout << "  Output file   : " << outRootPath << "\n";
  cout << "  Unfolding method        : " << m << "\n";
  if (m == "BAYES")
    cout << "  Bayes iters   : " << nIter << "\n\n";
  else
    cout << "  SVD k_reg     : ";
  if (m == "SVD") {
    for (int i = 0; i < kNSVDReg; ++i)
      cout << kSVDRegValues[i] << (i+1<kNSVDReg ? ", " : "\n\n");
  }

  // loop R, centrality, ptlead
  for (size_t iR = 0; iR < kRadii.size(); ++iR) {
    const string R = kRadii[iR];

    // choose area cut from R
    double Rval = 0.0;
    if (sscanf(R.c_str(), "R%lf", &Rval) != 1) {
      cout << "[warn] Could not parse radius from " << R << ", skipping.\n";
      continue;
    }
    double areaMin = 0.0;
    if      (Rval < 0.25) areaMin = CUT_AREA_02;
    else if (Rval < 0.35) areaMin = CUT_AREA_03;
    else                  areaMin = CUT_AREA_04;

    for (size_t iC = 0; iC < kCentralities.size(); ++iC) {
      const string C = kCentralities[iC];

      const string treePath = R + "/" + C + "/JetTree";
      TTree* tr = dynamic_cast<TTree*>(fData->Get(treePath.c_str()));
      if (!tr) {
        cout << "[note] missing tree in data: " << treePath << " (skip)\n";
        continue;
      }

      cout << "\n>>> Tag: " << R << "  " << C << "\n";

      // enable only needed branches (RECO only!)
      tr->SetBranchStatus("*", 0);
      tr->SetBranchStatus("reco_pt_corr", 1);
      tr->SetBranchStatus("reco_pt_lead", 1);
      tr->SetBranchStatus("reco_trigger_match", 1);
      tr->SetBranchStatus("centralityWeight", 1);
      tr->SetBranchStatus("reco_area", 1);
      tr->SetBranchStatus("reco_neutral_fraction", 1);

      float reco_pt_corr=0, reco_pt_lead=0;
      Bool_t reco_trigger_match=kFALSE;
      float centralityWeight=1.0f;
      float reco_area=0.0f, reco_neutral_fraction=0.0f;

      tr->SetBranchAddress("reco_pt_corr", &reco_pt_corr);
      tr->SetBranchAddress("reco_pt_lead", &reco_pt_lead);
      tr->SetBranchAddress("reco_trigger_match", &reco_trigger_match);
      tr->SetBranchAddress("centralityWeight", &centralityWeight);
      tr->SetBranchAddress("reco_area", &reco_area);
      tr->SetBranchAddress("reco_neutral_fraction", &reco_neutral_fraction);

      const Long64_t n = tr->GetEntries();
      if (n <= 0) {
        cout << "  [note] empty tree in data, skip.\n";
        continue;
      }

      for (int ic = 0; ic < kNPtLeadCuts; ++ic) {
        const double cut = kPtLeadCuts[ic];
        const string tag      = R + "_" + C + Form("_ptlead%.0f", cut);
        const string tagfile  = tag; // same as used when saving responses

        cout << "  - pTlead cut >= " << cut << " GeV (tag " << tag << ")\n";

        // ---------------- measured spectrum in DATA ----------------
        TH1D* hMeasData = new TH1D(
            ("hMeasData_"+tag).c_str(),
            ";p_{T}^{reco,corr} [GeV];counts",
            nbins_meas, bin_meas_edges);
        hMeasData->SetDirectory(0);

        // event loop
        for (Long64_t i = 0; i < n; ++i) {
          if ((i % 200000) == 0)
            cout << "    filling data: " << i << "/" << n << "\r" << std::flush;

          tr->GetEntry(i);

          const double w = (double)centralityWeight;

          if (!reco_trigger_match)               continue;
          if (reco_area < areaMin)               continue;
          if (reco_neutral_fraction > CUT_NEUTRAL_FRACTION) continue;
          if (reco_pt_lead < cut)                continue;

          hMeasData->Fill(reco_pt_corr, w);
        }
        cout << "\n    data integral = " << hMeasData->Integral(0,-1) << endl;

        // -------- Trigger efficiency correction (RECO binning, BEFORE unfolding) --------
         TH1D* hMeasData_raw = (TH1D*)hMeasData->Clone(Form("hMeasData_raw_%s", tag.c_str()));
         hMeasData_raw->SetDirectory(0);

        // // Load trigger efficiency (must match measured binning)
        // TH1D* hTrigEffReco = GetEffHistChecked(
        //     fEff,
        //     tagfile,
        //     "h_trig_eff_reco",
        //     hMeasData,
        //     Form("hTrigEffReco_%s", tag.c_str()),
        //     /*verbose=*/true);

        // if (hTrigEffReco) {
        //   ApplyDivideCorr(hMeasData, hTrigEffReco, /*verbose=*/true);
        //   cout << "    data integral after trig-eff corr = "
        //       << hMeasData->Integral(0,-1) << endl;
        // }

        TH1D* hTrigEffReco = 0; // TEMP: trig-eff disabled


        TH1D* hMeasData_trigCorr = (TH1D*)hMeasData->Clone(Form("hMeasData_trigCorr_%s", tag.c_str()));
        hMeasData_trigCorr->SetDirectory(0);

        // ---- Purity correction (fake removal) ----
        // TH1D* hPurityReco = GetEffHistChecked(
        //     fEff,
        //     tagfile,
        //     "h_pur_eff_reco",          
        //     hMeasData,                 // same binning as current measured spectrum
        //     Form("hPurityReco_%s", tag.c_str()),
        //     /*verbose=*/true);

        // if (hPurityReco) {
        //   ApplyMultiplyCorr(hMeasData, hPurityReco, /*verbose=*/true);
        //   cout << "    data integral after purity corr = "
        //       << hMeasData->Integral(0,-1) << endl;
        // }

         TH1D* hPurityReco = 0; // TEMP: purity correction disabled

        // ---------------- load response & prior ----------------
        TDirectory* dResp =
          dynamic_cast<TDirectory*>(fRespAll->Get(tagfile.c_str()));
        if (!dResp) {
          cout << "    [ERROR] no directory '" << tagfile << "' in response file --> skipping this tag\n";
          delete hMeasData;
          delete hMeasData_raw;
          delete hMeasData_trigCorr;
          if (hTrigEffReco) delete hTrigEffReco;
          if (hPurityReco)  delete hPurityReco;
          continue;
        }

        //   response      -> "response"
        //   prior (TH1D)  -> "hPrior"
        RooUnfoldResponse* response = 0;
        dResp->GetObject("response", response);
        if (!response) {
          cout << "    [ERROR] object 'response' not found in dir " << tagfile
              << " --> skipping\n";

          delete hMeasData;
          delete hMeasData_raw;
          delete hMeasData_trigCorr;

          if (hTrigEffReco) delete hTrigEffReco;
          if (hPurityReco)  delete hPurityReco;

          continue;
        }

        TH1D* hPrior = 0;
        dResp->GetObject("hPrior", hPrior);
        if (!hPrior) {
          cout << "    [WARN] 'hPrior' not found in dir " << tagfile
               << " (will use default prior from response)\n";
        }

       


        // =========================================================
        //                   BAYESIAN UNFOLDING
        // =========================================================

        if (m == "BAYES") {
        
          cout << "    running RooUnfoldBayes with " << nIter << " iterations\n";
          RooUnfoldBayes unfold(response, hMeasData, nIter,
                                /*smoothit=*/false,
                                hPrior); // may be 0, RooUnfold handles it

          TH1D* hUnfoldedRecoBins = (TH1D*)unfold.Hunfold(); // reco-binning
          hUnfoldedRecoBins->SetName(Form("UnfoldedRecoBins_%s", tag.c_str()));
          hUnfoldedRecoBins->SetDirectory(0);

          // Rebin to truth binning (this is what you'll use later)
          TH1D* hUnfoldedTruthBins_raw = (TH1D*)hUnfoldedRecoBins->Rebin(
              nbins_truth,
              Form("hUnfoldedTruthBins_raw_%s", tag.c_str()),
              bin_truth_edges);
          hUnfoldedTruthBins_raw->SetDirectory(0);

          // -------- Matching efficiency correction (TRUTH binning, AFTER unfolding) --------
          // TH1D* hMatchEffTruth = GetEffHistChecked(
          //     fEff,
          //     tagfile,
          //     "h_match_eff_truth",
          //     hUnfoldedTruthBins_raw,
          //     Form("hMatchEffTruth_%s", tag.c_str()),
          //     /*verbose=*/true);

          // TH1D* hUnfoldedTruthBins_matchCorr = 0;
          // if (hMatchEffTruth) {
          //   hUnfoldedTruthBins_matchCorr = (TH1D*)hUnfoldedTruthBins_raw->Clone(
          //       Form("hUnfoldedTruthBins_matchCorr_%s", tag.c_str()));
          //   hUnfoldedTruthBins_matchCorr->SetDirectory(0);

          //   ApplyDivideCorr(hUnfoldedTruthBins_matchCorr, hMatchEffTruth, /*verbose=*/true);
          // }

          TH1D* hMatchEffTruth = 0;               // TEMP: match-eff disabled
          TH1D* hUnfoldedTruthBins_matchCorr = 0; // keep null

          // ---------------- save to output file ----------------
          TDirectory* dOut = fOutAll->mkdir(tagfile.c_str());
          if (!dOut) dOut = fOutAll->GetDirectory(tagfile.c_str());
          dOut->cd();

          hMeasData_raw->Write("hMeasData_raw", TObject::kOverwrite);
          hMeasData_trigCorr->Write("hMeasData_trigCorr", TObject::kOverwrite);
          hMeasData->Write("hMeasData_trigPurCorr", TObject::kOverwrite);
          if (hTrigEffReco)  hTrigEffReco->Write("hTrigEffReco", TObject::kOverwrite);
          if (hPurityReco)   hPurityReco->Write("hPurityReco", TObject::kOverwrite);
          hUnfoldedRecoBins->Write("hUnfoldedRecoBins", TObject::kOverwrite);
          hUnfoldedTruthBins_raw->Write("hUnfoldedTruthBins_raw", TObject::kOverwrite);
          if (hMatchEffTruth) hMatchEffTruth->Write("hMatchEffTruth", TObject::kOverwrite);
          if (hUnfoldedTruthBins_matchCorr) hUnfoldedTruthBins_matchCorr->Write("hUnfoldedTruthBins_matchCorr", TObject::kOverwrite);

          TMatrixD cov = unfold.Eunfold(RooUnfold::kCovariance);
          cov.Write("covariance", TObject::kOverwrite);
          
          // cleanup
          delete hMeasData;
          delete hMeasData_raw;
          if (hTrigEffReco) delete hTrigEffReco;

          delete hUnfoldedRecoBins;
          delete hUnfoldedTruthBins_raw;

          delete hMeasData_trigCorr;
          if (hPurityReco) delete hPurityReco;
          if (hMatchEffTruth) delete hMatchEffTruth;
          if (hUnfoldedTruthBins_matchCorr) delete hUnfoldedTruthBins_matchCorr;
        } //--- end of BAYESIAN unfolding ---


        // =========================================================
        //                     SVD UNFOLDING
        // =========================================================
        else if (m == "SVD") {
          cout << "    running RooUnfoldSvd for k_reg = ";
          for (int ir = 0; ir < kNSVDReg; ++ir)
            cout << kSVDRegValues[ir] << (ir+1<kNSVDReg ? ", " : "\n");
 
          vector<TH1D*> unfoldedSVD(kNSVDReg, nullptr);
          vector<TH1D*> unfoldedTruthSVD(kNSVDReg, nullptr);
 
          for (int ir = 0; ir < kNSVDReg; ++ir) {
            const int reg = kSVDRegValues[ir];
 
            RooUnfoldSvd unfoldSVD(response, hMeasData, reg);
            TH1D* hUnfoldedRecoBinsSVD = (TH1D*)unfoldSVD.Hunfold();
            if (!hUnfoldedRecoBinsSVD) {
              cout << "    [ERROR] SVD unfolding failed for reg=" << reg << "\n";
              continue;
            }
            hUnfoldedRecoBinsSVD->SetDirectory(nullptr);
            hUnfoldedRecoBinsSVD->SetName(Form("UnfoldedSVD_%s_reg%d", tag.c_str(), reg));
            unfoldedSVD[ir] = hUnfoldedRecoBinsSVD;
 
            // --- print diagnostics (same as embedding code) ---
            auto* svdImpl = unfoldSVD.Impl();
            if (svdImpl) {
              TVectorD sv = svdImpl->GetSV();
              cout << "=== Singular values (reg=" << reg << ") ===" << endl;
              for (int i = 0; i < sv.GetNrows(); i++)
                cout << Form("  Mode %2d : SV = %10.6e", i, sv[i]) << endl;
 
              TH1* d = svdImpl->GetD();
              if (d) {
                cout << "=== d-vector (reg=" << reg << ") ===" << endl;
                for (int i = 1; i <= d->GetNbinsX(); i++)
                  cout << Form("  Mode %2d : |d_i| = %10.6e", i, fabs(d->GetBinContent(i))) << endl;
              } else {
                cout << "    [note] d still null after Hunfold() for reg=" << reg << "\n";
              }
            }
            // rebin to truth binning
            unfoldedTruthSVD[ir] = (TH1D*)hUnfoldedRecoBinsSVD->Rebin(
                nbins_truth,
                Form("UnfSVDTruthBins_%d_%s", reg, tag.c_str()),
                bin_truth_edges);
            if (unfoldedTruthSVD[ir]) unfoldedTruthSVD[ir]->SetDirectory(nullptr);
          }
                  //  -------- Matching efficiency correction (TRUTH binning, AFTER unfolding) --------
          TH1D* hMatchEffTruth = GetEffHistChecked(
              fEff,
              tagfile,
              "h_match_eff_truth",
              unfoldedTruthSVD[0], // use first SVD reg as reference for binning
              Form("hMatchEffTruth_%s", tag.c_str()),
              /*verbose=*/true);

          TH1D* hUnfoldedTruthBins_matchCorr = 0;
          if (hMatchEffTruth) {
            hUnfoldedTruthBins_matchCorr = (TH1D*)unfoldedTruthSVD[0]->Clone(
                Form("hUnfoldedTruthBins_matchCorr_%s", tag.c_str()));
            hUnfoldedTruthBins_matchCorr->SetDirectory(0);

            ApplyDivideCorr(hUnfoldedTruthBins_matchCorr, hMatchEffTruth, /*verbose=*/true);
          }

          //TH1D* hMatchEffTruth = 0;               // TEMP: match-eff disabled
          //TH1D* hUnfoldedTruthBins_matchCorr = 0; // keep null

          // ---------------- save to output file ----------------
          TDirectory* dOut = fOutAll->mkdir(tagfile.c_str());
          if (!dOut) dOut = fOutAll->GetDirectory(tagfile.c_str());
          dOut->cd();

          hMeasData_raw->Write("hMeasData_raw", TObject::kOverwrite);
          hMeasData_trigCorr->Write("hMeasData_trigCorr", TObject::kOverwrite);
          hMeasData->Write("hMeasData_trigPurCorr", TObject::kOverwrite);
          if (hTrigEffReco)  hTrigEffReco->Write("hTrigEffReco", TObject::kOverwrite);
          if (hPurityReco)   hPurityReco->Write("hPurityReco", TObject::kOverwrite);
            for (int ir = 0; ir < kNSVDReg; ++ir) {
            const int reg = kSVDRegValues[ir];
            if (unfoldedSVD[ir])
              unfoldedSVD[ir]->Write(Form("hUnfoldedRecoBinsSVD_reg%d", reg), TObject::kOverwrite);
            if (unfoldedTruthSVD[ir])
              unfoldedTruthSVD[ir]->Write(Form("hUnfoldedTruthBins_rawSVD_reg%d", reg), TObject::kOverwrite);
          }
          if (hMatchEffTruth) hMatchEffTruth->Write("hMatchEffTruth", TObject::kOverwrite);
          if (hUnfoldedTruthBins_matchCorr) hUnfoldedTruthBins_matchCorr->Write("hUnfoldedTruthBins_matchCorr", TObject::kOverwrite);

          // covariance for the first regularization value
          if (unfoldedSVD[0]) {
            const int regDefault = kSVDRegValues[0];
            RooUnfoldSvd uCov(response, hMeasData, regDefault);
            uCov.Hunfold();
            TMatrixD cov = uCov.Eunfold(RooUnfold::kCovariance);
            cov.Write(Form("covariance_reg%d", regDefault), TObject::kOverwrite);
          }

          // cleanup
          delete hMeasData;
          delete hMeasData_raw;
          if (hTrigEffReco) delete hTrigEffReco;

          for (int ir = 0; ir < kNSVDReg; ++ir) {
            if (unfoldedSVD[ir])      delete unfoldedSVD[ir];
            if (unfoldedTruthSVD[ir]) delete unfoldedTruthSVD[ir];
          }

          delete hMeasData_trigCorr;
          if (hPurityReco) delete hPurityReco;
          if (hMatchEffTruth) delete hMatchEffTruth;
          if (hUnfoldedTruthBins_matchCorr) delete hUnfoldedTruthBins_matchCorr;
        } //--- end of SVD unfolding ---
      
 


      } // ptlead cuts
    } // centralities
  } // radii

  fOutAll->Write();
  fOutAll->Close();
  delete fOutAll;

  fRespAll->Close();
  delete fRespAll;

  fData->Close();
  delete fData;

  fEff->Close();
  delete fEff;

  cout << "\nAll done. Unfolded data written to: "
       << outRootPath << "\n";
}

// ==================================================
// Helper implementations
// ==================================================

static bool SameBinning1D(const TH1* a, const TH1* b, double tol)
{
  if (!a || !b) return false;
  if (a->GetNbinsX() != b->GetNbinsX()) return false;

  const TAxis* ax = a->GetXaxis();
  const TAxis* bx = b->GetXaxis();

  for (int ib = 1; ib <= a->GetNbinsX(); ++ib) {
    if (std::fabs(ax->GetBinLowEdge(ib) - bx->GetBinLowEdge(ib)) > tol) return false;
    if (std::fabs(ax->GetBinUpEdge(ib)  - bx->GetBinUpEdge(ib))  > tol) return false;
  }
  return true;
}

static TH1D* GetEffHistChecked(TFile* fEff,
                               const std::string& tagdirName,
                               const char* objName,
                               const TH1* ref,
                               const char* cloneNameForOut,
                               bool verbose)
{
  if (!fEff || !ref) return 0;

  TDirectory* d = dynamic_cast<TDirectory*>(fEff->Get(tagdirName.c_str()));
  if (!d) {
    if (verbose) cout << "    [ERROR] no dir '" << tagdirName << "' in eff file\n";
    return 0;
  }

  TH1D* h = 0;
  d->GetObject(objName, h);
  if (!h) {
    if (verbose) cout << "    [ERROR] missing '" << objName << "' in " << tagdirName << "\n";
    return 0;
  }

  if (!SameBinning1D(h, ref)) {
    if (verbose) cout << "    [ERROR] binning mismatch for " << objName << "\n";
    return 0;
  }

  TH1D* hc = (TH1D*)h->Clone(cloneNameForOut);
  hc->SetDirectory(0);
  return hc;
}

static void ApplyDivideCorr(TH1D* spec,
                            const TH1D* eff,
                            bool verbose)
{
  if (!spec || !eff) return;

  for (int ib = 1; ib <= spec->GetNbinsX(); ++ib) {
    const double y  = spec->GetBinContent(ib);
    const double ey = spec->GetBinError(ib);
    const double e  = eff->GetBinContent(ib);
    const double ee = eff->GetBinError(ib);

    if (e > 0.0) {
      const double ycorr = y / e;
      double ecorr = (y > 0.0)
        ? ycorr * std::sqrt((ey*ey)/(y*y) + (ee*ee)/(e*e))
        : ey / e;

      spec->SetBinContent(ib, ycorr);
      spec->SetBinError(ib, ecorr);
    } else {
      spec->SetBinContent(ib, 0.0);
      spec->SetBinError(ib, 0.0);
    }
  }

  if (verbose)
    cout << "    applied divide correction: " << eff->GetName() << endl;
}

static void ApplyMultiplyCorr(TH1D* spec,
                              const TH1D* fac,
                              bool verbose)
{
  if (!spec || !fac) return;

  for (int ib = 1; ib <= spec->GetNbinsX(); ++ib) {
    const double y  = spec->GetBinContent(ib);
    const double ey = spec->GetBinError(ib);
    const double f  = fac->GetBinContent(ib);
    const double ef = fac->GetBinError(ib);

    // Multiplication: y' = y * f
    const double ycorr = y * f;

    double ecorr = 0.0;
    if (y != 0.0 && f != 0.0) {
      const double rel2 = (ey*ey)/(y*y) + (ef*ef)/(f*f);
      ecorr = std::fabs(ycorr) * std::sqrt(rel2);
    } else {
      // if one factor is zero, absolute propagation:
      // ycorr = y*f, so d(ycorr)^2 = (f*ey)^2 + (y*ef)^2
      ecorr = std::sqrt((f*ey)*(f*ey) + (y*ef)*(y*ef));
    }

    spec->SetBinContent(ib, ycorr);
    spec->SetBinError(ib, ecorr);
  }

  if (verbose)
    cout << "    applied multiply correction: " << fac->GetName() << endl;
}