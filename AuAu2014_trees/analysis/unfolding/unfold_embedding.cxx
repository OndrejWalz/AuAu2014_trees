#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"

#include <vector>
#include <string>
#include <iostream>

using std::string;
using std::vector;
using std::cout;
using std::endl;

// ------------------------- config -----------------------------

static const int kBayesIters[] = {1, 4, 5, 7};
static const int kNBayesIters = sizeof(kBayesIters)/sizeof(kBayesIters[0]);

static const double kTestFrac = 0.50;     // 50/50 split
static const UInt_t kSeed     = 12345;    // deterministic split

static const double kPtLeadCuts[] = {0.0, 5.0, 7.0, 9.0};
static const int    kNPtLeadCuts  = sizeof(kPtLeadCuts)/sizeof(kPtLeadCuts[0]);

// measured & truth binning
static const int nbins_meas = 24;
static const double bin_meas_edges[nbins_meas+1] = {
  -100,-80,-60,-40,-20,-10,-5,-2.5,0,2.5,5,7.5,10,12.5,15,17.5,
  20,22.5,25,27.5,30,35,40,50,60
};

static const int nbins_truth = 10;
static const double bin_truth_edges[nbins_truth+1] = {
  0,5,10,15,20,25,30,35,40,50,60
};

static const vector<string> kCentralities =
  {"CENT_0_10", "MID_20_40", "PERI_60_80"};
static const vector<string> kRadii =
  {"R0.2", "R0.3", "R0.4"};

// ---- pThat bins (upper edges) and xsec weights (same order) ----
static const int kNPthatBins = 11;

// upper edges for each bin: 3_5 -> 5, ..., 40_50 -> 50, 50_-1 -> -1 (no upper bound)
static const double kPtHatMax[kNPthatBins] =
  {5, 7, 9, 11, 15, 20, 25, 30, 40, 50, -1};

// xsecWeight values that tag each bin (from your code)
static const double kXsecWeights[kNPthatBins] =
  {1.616e+0,  1.355e-01, 2.288e-02, 5.524e-03, 2.203e-03,
   3.437e-04, 4.681e-05, 8.532e-06, 2.178e-06, 1.198e-07, 6.939e-09};

// reco dummy sentinel (keep real negative jets, reject dummy ~ -999)
static const double RECO_PTCORR_DUMMY_CUT = -500.0;

// match xsecWeight -> bin index (tolerant compare)
static int FindPtHatBin(double xsecW)
{
  // relative tolerance, because floats are a joy
  const double relTol = 1e-6;
  for (int i = 0; i < kNPthatBins; ++i) {
    double ref = kXsecWeights[i];
    double diff = fabs(xsecW - ref);
    if (diff <= relTol * fabs(ref)) return i;
  }
  return -1; // unknown
}

// ------------------ jet-quality cuts -------------------------

// Area cuts per jet radius
static const double CUT_AREA_02 = 0.07;  // R = 0.2
static const double CUT_AREA_03 = 0.20;  // R = 0.3
static const double CUT_AREA_04 = 0.40;  // R = 0.4

// Max neutral energy fraction
static const double CUT_NEUTRAL_FRACTION = 0.95;

// --------------------------------------------------------------

static void EnsureDir(const string& path){
  if (gSystem->AccessPathName(path.c_str()))
    gSystem->mkdir(path.c_str(), /*recursive=*/true);
}

// Main macro: builds closure + full-statistics responses
void unfold_embedding(const char* inputFile,
                      const char* outDir)
{
  gStyle->SetOptStat(0);
  EnsureDir(outDir);
  TH1::SetDefaultSumw2(kTRUE);

  TFile* fin = TFile::Open(inputFile, "READ");
  if (!fin || fin->IsZombie()) {
    cout << "[error] Cannot open input file: " << inputFile << endl;
    return;
  }

  const string outRootPath = string(outDir) + "/responses_embedding.root";
  TFile* fout = TFile::Open(outRootPath.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    cout << "[error] Cannot create output file: " << outRootPath << endl;
    fin->Close(); delete fin;
    return;
  }

  cout << "Input embedding : " << inputFile  << endl;
  cout << "Output responses: " << outRootPath << endl;

  // loop R, centrality, ptlead
  for (size_t iR = 0; iR < kRadii.size(); ++iR) {
    const string R = kRadii[iR];

    // parse numeric R to choose area cut
    double Rval = 0.0;
    if (sscanf(R.c_str(), "R%lf", &Rval) != 1) {
      cout << "[warn] Could not parse radius from " << R << ", skipping.\n";
      continue;
    }
    double areaMin = 0.0;
    if      (Rval < 0.25) areaMin = CUT_AREA_02; // ~0.2
    else if (Rval < 0.35) areaMin = CUT_AREA_03; // ~0.3
    else                  areaMin = CUT_AREA_04; // ~0.4

    for (size_t iC = 0; iC < kCentralities.size(); ++iC) {
      const string C = kCentralities[iC];

      const string treePath = R + "/" + C + "/JetTree";
      TTree* tr = dynamic_cast<TTree*>(fin->Get(treePath.c_str()));
      if (!tr) {
        cout << "[note] missing tree: " << treePath << " (skip)\n";
        continue;
      }

      // speed up I/O: enable only needed branches
      tr->SetBranchStatus("*", 0);
      tr->SetBranchStatus("mc_pt", 1);
      tr->SetBranchStatus("mc_pt_lead", 1);
      tr->SetBranchStatus("reco_pt_corr", 1);
      tr->SetBranchStatus("reco_pt_lead", 1);
      tr->SetBranchStatus("reco_trigger_match", 1);
      tr->SetBranchStatus("centralityWeight", 1);
      tr->SetBranchStatus("xsecWeight", 1);
      tr->SetBranchStatus("reco_area", 1);
      tr->SetBranchStatus("reco_neutral_fraction", 1);

      // branch addresses
      float mc_pt=0, mc_pt_lead=0;
      float reco_pt_corr=0, reco_pt_lead=0;
      bool  reco_trigger_match=false;
      float centralityWeight=1.0f, xsecWeight=1.0f;
      float reco_area=0.0f, reco_neutral_fraction=0.0f;

      tr->SetBranchAddress("mc_pt", &mc_pt);
      tr->SetBranchAddress("mc_pt_lead", &mc_pt_lead);
      tr->SetBranchAddress("reco_pt_corr", &reco_pt_corr);
      tr->SetBranchAddress("reco_pt_lead", &reco_pt_lead);
      tr->SetBranchAddress("reco_trigger_match", &reco_trigger_match);
      tr->SetBranchAddress("centralityWeight", &centralityWeight);
      tr->SetBranchAddress("xsecWeight", &xsecWeight);
      tr->SetBranchAddress("reco_area", &reco_area);
      tr->SetBranchAddress("reco_neutral_fraction", &reco_neutral_fraction);

      const Long64_t n = tr->GetEntries();
      if (n <= 0) {
        cout << "[note] empty tree: " << treePath << " (skip)\n";
        continue;
      }

      cout << "\n=== " << R << "  " << C << " (entries: " << n << ") ===\n";

      for (int ic = 0; ic < kNPtLeadCuts; ++ic) {
        const double cut = kPtLeadCuts[ic];
        const string tag     = R + "_" + C + Form("_ptlead%.0f", cut);
        const string tagfile = tag;

        cout << "  >> pTlead >= " << cut << " GeV  (tag " << tag << ")\n";

        // --- histos: train/test for closure ---
        TH1D* hMeasTrain = new TH1D(("hMeasTrain_"+tag).c_str(),
            ";reco p_{T}^{corr} [GeV];dN/dp_{T}",
            nbins_meas, bin_meas_edges);
        TH1D* hTrueTrain = new TH1D(("hTrueTrain_"+tag).c_str(),
            ";mc p_{T} [GeV];dN/dp_{T}",
            nbins_truth, bin_truth_edges);
        TH1D* hMeasTest  = (TH1D*)hMeasTrain->Clone(("hMeasTest_"+tag).c_str());
        TH1D* hTrueTest  = (TH1D*)hTrueTrain->Clone(("hTrueTest_"+tag).c_str());

        hMeasTrain->SetDirectory(0);
        hTrueTrain->SetDirectory(0);
        hMeasTest ->SetDirectory(0);
        hTrueTest ->SetDirectory(0);

        // --- full-stat histos (no split) ---
        TH1D* hMeasFull = new TH1D(("hMeasFull_"+tag).c_str(),
            ";reco p_{T}^{corr} [GeV];dN/dp_{T}",
            nbins_meas, bin_meas_edges);
        TH1D* hTrueFull = new TH1D(("hTrueFull_"+tag).c_str(),
            ";mc p_{T} [GeV];dN/dp_{T}",
            nbins_truth, bin_truth_edges);
        hMeasFull->SetDirectory(0);
        hTrueFull->SetDirectory(0);

        // --- response matrices ---
        // For closure (train only)
        TH2D* hRespTrain = new TH2D(("hRespTrain_"+tag).c_str(),
            ";p_{T}^{reco,corr};p_{T}^{mc}",
            nbins_meas, bin_meas_edges,
            nbins_truth, bin_truth_edges);
        hRespTrain->SetDirectory(0);

        // For full-statistics response (train+test together)
        TH2D* hRespFull = new TH2D(("hRespFull_"+tag).c_str(),
            ";p_{T}^{reco,corr};p_{T}^{mc}",
            nbins_meas, bin_meas_edges,
            nbins_truth, bin_truth_edges);
        hRespFull->SetDirectory(0);

        // --- prior (truth-only, using *all* jets passing truth-side cuts) ---
        TH1D* hPrior = new TH1D(("hPrior_"+tag).c_str(),
            ";mc p_{T} [GeV];prior",
            nbins_truth, bin_truth_edges);
        hPrior->SetDirectory(0);

        TRandom3 rng(kSeed);

        // event loop
        for (Long64_t i = 0; i < n; ++i) {
          if ((i % 200000) == 0)
            cout << "    ["<< tag <<"] " << i << "/" << n << "\r" << std::flush;

          tr->GetEntry(i);

          // weight (same for prior and response)
          const double w = (double)centralityWeight * (double)xsecWeight;

          const bool haveMC = (mc_pt > 0.0);
          const bool haveReco = (reco_pt_corr > RECO_PTCORR_DUMMY_CUT);

          // ----- fill prior: only MC-side cuts -----
          if (haveMC && mc_pt_lead >= cut) {
            hPrior->Fill(mc_pt, w);
          }

          // ----- reco-side cuts for response & closure -----
          if (!haveReco) continue;

          // existing quality/trigger cuts
          if (!reco_trigger_match) continue;
          if (reco_area < areaMin) continue;

          // dual ptlead cut (both reco & MC)
          if (!(reco_pt_lead >= cut && mc_pt_lead >= cut)) continue;

          // ===== full-statistics response (no split) =====
          hRespFull->Fill(reco_pt_corr, mc_pt, w);
          hMeasFull->Fill(reco_pt_corr, w);
          hTrueFull->Fill(mc_pt,        w);

          // ===== closure split: train vs test =====
          const bool train = (rng.Uniform() > kTestFrac);
          if (train) {
            hRespTrain->Fill(reco_pt_corr, mc_pt, w);
            hMeasTrain->Fill(reco_pt_corr, w);
            hTrueTrain->Fill(mc_pt,        w);
          } else {
            hMeasTest->Fill(reco_pt_corr, w);
            hTrueTest->Fill(mc_pt,        w);
          }
        }
        cout << endl;

        // normalize prior to training truth (as before)
        double intPrior = hPrior->Integral();
        double intTrue  = hTrueTrain->Integral();
        if (intPrior > 0.0 && intTrue > 0.0) {
          hPrior->Scale(intTrue / intPrior);
        }

        cout << "    Integral(Truth train) = " << hTrueTrain->Integral(0, -1) << endl;
        cout << "    Integral(Meas  train) = " << hMeasTrain->Integral(0, -1) << endl;
        cout << "    Integral(Truth test ) = " << hTrueTest ->Integral(0, -1) << endl;
        cout << "    Integral(Meas  test ) = " << hMeasTest ->Integral(0, -1) << endl;
        cout << "    Integral(Prior      ) = " << hPrior    ->Integral(0, -1) << endl;
        cout << "    Integral(Truth full)  = " << hTrueFull ->Integral(0, -1) << endl;
        cout << "    Integral(Meas  full)  = " << hMeasFull ->Integral(0, -1) << endl;

        // --- build responses ---

        // Full-statistics response: this one you will use for unfolding DATA
        RooUnfoldResponse response_full(hMeasFull, hTrueFull, hRespFull);
        response_full.SetName(("response_full_"+tag).c_str());

        // Closure response (train-only): used only for unfolding test sample
        RooUnfoldResponse response_closure(hMeasTrain, hTrueTrain, hRespTrain);
        response_closure.SetName(("response_closure_"+tag).c_str());

        // --- unfolding (closure) with explicit prior ---
        vector<TH1D*> unfolded(kNBayesIters, 0);
        for (int ib = 0; ib < kNBayesIters; ++ib) {
          RooUnfoldBayes u(&response_closure, hMeasTest,
                           kBayesIters[ib], false, hPrior);
          TH1D* hunf = (TH1D*)u.Hunfold();
          hunf->SetDirectory(0);
          hunf->SetName(Form("Unfolded_%s_iter%d", tag.c_str(), kBayesIters[ib]));
          unfolded[ib] = hunf;
        }

        // Rebin unfolded to truth binning (for plotting & ratio)
        vector<TH1D*> unfoldedTruth(kNBayesIters, 0);
        for (int ib = 0; ib < kNBayesIters; ++ib) {
          if (!unfolded[ib]) continue;
          unfoldedTruth[ib] = (TH1D*)unfolded[ib]->Rebin(
              nbins_truth,
              Form("UnfTruthBins_%d_%s", kBayesIters[ib], tag.c_str()),
              bin_truth_edges);
          unfoldedTruth[ib]->SetDirectory(0);
        }

        // ========================= PLOTTING (closure) =========================
        TCanvas* c = new TCanvas(("c_"+tag).c_str(), "", 800, 1200);
        c->Divide(1,2);

        // ---------- top pad: shapes in TRUTH binning ----------
        c->cd(1);
        gPad->SetLogy();

        TH1D* hT_plot = (TH1D*)hTrueTest->Clone(("hTrueW_"+tag).c_str());
        hT_plot->SetMarkerStyle(20);
        hT_plot->SetLineColor(kBlack);

        TH1D* hM_truth = (TH1D*)hMeasTest->Rebin(
            nbins_truth,
            ("hMeasTruthBins_"+tag).c_str(),
            bin_truth_edges);
        hM_truth->SetMarkerStyle(24);
        hM_truth->SetLineColor(kBlue+1);

        hT_plot->GetXaxis()->SetTitle("p_{T} [GeV]");
        hT_plot->Draw("E1");
        hM_truth->Draw("E1 SAME");

        TLegend* leg = new TLegend(0.58,0.58,0.88,0.88);
        leg->AddEntry(hT_plot,  "Truth (test)",    "lp");
        leg->AddEntry(hM_truth, "Measured (test)", "lp");

        for (int ib = 0; ib < kNBayesIters; ++ib) {
          TH1D* w = unfoldedTruth[ib];
          if (!w) continue;
          w->SetMarkerStyle(20);
          w->SetMarkerColor(kMagenta + ib);
          w->SetLineColor(kMagenta + ib);
          w->Draw("E1 SAME");
          leg->AddEntry(w, Form("Bayes %d it.", kBayesIters[ib]), "lp");
        }
        leg->Draw();

        // ---------- bottom pad: ratio unfolded / truth ----------
        c->cd(2);
        gPad->SetGridy();

        TH1D* firstRatio = 0;

        for (int ib = 0; ib < kNBayesIters; ++ib) {
          TH1D* hunf_reb = unfoldedTruth[ib];
          if (!hunf_reb) continue;

          TH1D* r = (TH1D*)hunf_reb->Clone(
              Form("ratio_%d_%s", kBayesIters[ib], tag.c_str()));
          r->SetDirectory(0);
          r->Divide(hTrueTest);

          if (!firstRatio) {
            firstRatio = r;
            firstRatio->SetTitle("");
            firstRatio->GetYaxis()->SetTitle("Unfolded / Truth");
            firstRatio->GetYaxis()->SetRangeUser(0.0, 2.0);
            firstRatio->Draw("E1");
          } else {
            r->Draw("E1 SAME");
          }
        }

        if (firstRatio) {
          double xmin = firstRatio->GetXaxis()->GetXmin();
          double xmax = firstRatio->GetXaxis()->GetXmax();

          TLine* ln1 = new TLine(xmin, 1.1, xmax, 1.1);
          ln1->SetLineStyle(2); ln1->SetLineColor(kGray+1); ln1->Draw();

          TLine* ln2 = new TLine(xmin, 0.9, xmax, 0.9);
          ln2->SetLineStyle(2); ln2->SetLineColor(kGray+1); ln2->Draw();
        }

        // ---------- save ----------
        const string pdfPath = string(outDir) + "/closure_" + tagfile + ".pdf";

        // create a directory for this (R, centrality, ptlead) inside fout
        TDirectory* d = fout->mkdir(tagfile.c_str());
        if (!d) d = fout->GetDirectory(tagfile.c_str());
        d->cd();

        // Store full-statistics response (for data)
        hRespFull->Write("hRespRecoVsTruth_full");
        response_full.Write("response");  // <-- THIS is what unfold_data reads

        hMeasFull->Write("hMeasFull");
        hTrueFull->Write("hTrueFull");

        // Store closure response + components
        hRespTrain->Write("hRespRecoVsTruth_train");
        response_closure.Write("response_closure");

        hMeasTrain->Write("hMeasTrain");
        hTrueTrain->Write("hTrueTrain");
        hMeasTest ->Write("hMeasTest");
        hTrueTest ->Write("hTrueTest");
        hPrior    ->Write("hPrior");  // save the prior used

        // optionally save unfolded spectra for closure checks
        for (int ib = 0; ib < kNBayesIters; ++ib) {
          if (unfolded[ib]) {
            unfolded[ib]->Write(
              Form("Unfolded_iter%d", kBayesIters[ib])
            );
          }
        }

        c->SaveAs(pdfPath.c_str());

        // ---------- cleanup ----------
        delete c;
        delete hT_plot;
        delete hM_truth;
        delete hRespTrain;
        delete hRespFull;
        delete hMeasTrain;
        delete hTrueTrain;
        delete hMeasTest;
        delete hTrueTest;
        delete hMeasFull;
        delete hTrueFull;
        delete hPrior;
        delete leg;
        for (size_t ib = 0; ib < unfolded.size(); ++ib) {
          if (unfolded[ib]) delete unfolded[ib];
        }
        for (size_t ib = 0; ib < unfoldedTruth.size(); ++ib) {
          if (unfoldedTruth[ib]) delete unfoldedTruth[ib];
        }
      } // ptlead cuts
    } // centralities
  } // radii

  fout->Write();
  fout->Close();
  delete fout;

  fin->Close();
  delete fin;

  cout << "\nDone. Responses (closure + full) written to: "
       << outRootPath << endl;
}
