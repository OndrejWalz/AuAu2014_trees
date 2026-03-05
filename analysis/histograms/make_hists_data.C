#include "TFile.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TKey.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"

#include <iostream>
#include <string>
#include <cmath>

// =========================================================
//  Jet-quality cuts (must match maker!)
// =========================================================
const double CUT_AREA_02 = 0.07;  // R = 0.2
const double CUT_AREA_03 = 0.20;  // R = 0.3
const double CUT_AREA_04 = 0.40;  // R = 0.4
const double CUT_NEUTRAL_FRACTION = 0.95;
const double RECO_DUMMY_CUT = -500.0;  // your dummy threshold

// pTlead thresholds
const int    N_LEAD = 4;
const double PTLEAD_THR[N_LEAD] = {0.0, 5.0, 7.0, 9.0};

// helper: choose area cut based on R
double areaMinForR(double R) {
  if (R < 0.25) return CUT_AREA_02;
  if (R < 0.35) return CUT_AREA_03;
  return CUT_AREA_04;
}

void make_hists_data(const char *infile  = "data_merged.root",
                     const char *outfile = "hists_data.root")
{
  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetOptStat(0);
  gROOT->SetBatch(kTRUE);

  // where to write PDFs
  const char* pdfTop = "pdf";
  gSystem->mkdir(pdfTop, kTRUE);

  TFile *fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Error: cannot open input file " << infile << std::endl;
    return;
  }

  TFile *fout = TFile::Open(outfile, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "Error: cannot create output file " << outfile << std::endl;
    fin->Close();
    return;
  }

  std::cout << "Input : " << infile  << std::endl;
  std::cout << "Output: " << outfile << std::endl;

  // Loop over R directories: R0.2, R0.3, R0.4, ...
  TIter nextR(fin->GetListOfKeys());
  TKey *keyR = 0;

  while ((keyR = (TKey*) nextR())) {
    if (std::string(keyR->GetClassName()) != "TDirectoryFile") continue;

    std::string rname = keyR->GetName();
    if (rname.empty() || rname[0] != 'R') continue;

    double R = 0.0;
    if (sscanf(rname.c_str(), "R%lf", &R) != 1) continue;

    const double areaMin = areaMinForR(R);

    TDirectoryFile *dirR = (TDirectoryFile*) keyR->ReadObj();
    if (!dirR) continue;

    std::cout << "Radius directory: " << rname
              << " (R=" << R << ", areaMin=" << areaMin << ")\n";

    // output dir R
    fout->cd();
    TDirectory *outR = fout->mkdir(rname.c_str());
    if (!outR) outR = fout->GetDirectory(rname.c_str());

    // pdf dir for this R
    TString pdfRDir; pdfRDir.Form("%s/%s", pdfTop, rname.c_str());
    gSystem->mkdir(pdfRDir.Data(), kTRUE);

    // Loop over centrality subdirs
    TIter nextC(dirR->GetListOfKeys());
    TKey *keyC = 0;

    while ((keyC = (TKey*) nextC())) {
      if (std::string(keyC->GetClassName()) != "TDirectoryFile") continue;

      std::string cname = keyC->GetName();
      TDirectoryFile *dirC = (TDirectoryFile*) keyC->ReadObj();
      if (!dirC) continue;

      TTree *tree = (TTree*) dirC->Get("JetTree");
      if (!tree) continue;

      Long64_t nentries = tree->GetEntries();
      std::cout << "  Centrality: " << cname << " (entries: " << nentries << ")\n";

      // Require data-like branches (no MC/xsec/deltaR requirements)
      if (!tree->GetBranch("reco_pt") ||
          !tree->GetBranch("reco_pt_corr") ||
          !tree->GetBranch("reco_area") ||
          !tree->GetBranch("reco_neutral_fraction") ||
          !tree->GetBranch("reco_trigger_match") ||
          !tree->GetBranch("reco_pt_lead")) {
        std::cerr << "    Missing required reco branches in " << rname << "/" << cname
                  << " -> skipping.\n";
        continue;
      }

      // Optional weight branch (present in your embedding trees; may or may not exist in data trees)
        Bool_t HAS_CENT_WEIGHT = kFALSE;
        if (tree->GetBranch("centralityWeight") != 0) HAS_CENT_WEIGHT = kTRUE;
      // output dirs
      outR->cd();
      TDirectory *outC = outR->mkdir(cname.c_str());
      if (!outC) outC = outR->GetDirectory(cname.c_str());
      outC->cd();

      // pdf dir for this centrality
      TString pdfCDir; pdfCDir.Form("%s/%s", pdfRDir.Data(), cname.c_str());
      gSystem->mkdir(pdfCDir.Data(), kTRUE);

      // ---------------- branches ----------------
      Float_t centralityWeight = 1.f;

      Float_t reco_pt = -999.f;
      Float_t reco_pt_corr = -999.f;
      Float_t reco_area = -999.f;
      Float_t reco_neutral_fraction = 999.f;
      Float_t reco_pt_lead = -999.f;
      Bool_t  reco_trigger_match = kFALSE;

      if (HAS_CENT_WEIGHT) {
        tree->SetBranchAddress("centralityWeight", &centralityWeight);
      }

      tree->SetBranchAddress("reco_pt",                &reco_pt);
      tree->SetBranchAddress("reco_pt_corr",           &reco_pt_corr);
      tree->SetBranchAddress("reco_area",              &reco_area);
      tree->SetBranchAddress("reco_neutral_fraction",  &reco_neutral_fraction);
      tree->SetBranchAddress("reco_pt_lead",           &reco_pt_lead);
      tree->SetBranchAddress("reco_trigger_match",     &reco_trigger_match);

      // ---------------- histogram booking ----------------
      const int nbins_corr = 160;
      const double xmin_corr = -40.0;
      const double xmax_corr =  60.0;

      // spectra per ptlead threshold (reco_pt_corr)
      TH1D* hSpec[N_LEAD];
      for (int it = 0; it < N_LEAD; ++it) {
        hSpec[it] = new TH1D(
          Form("hSpec_recoPtCorr_ptlead%.0f_%s_%s", PTLEAD_THR[it], rname.c_str(), cname.c_str()),
          Form("Reco p_{T}^{corr} spectrum (ptlead>=%.0f);p_{T}^{reco,corr} [GeV];counts",
               PTLEAD_THR[it]),
          nbins_corr, xmin_corr, xmax_corr
        );
      }

      // ---------------- loop entries ----------------
      for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        // weight (data: no xsecWeight; use centralityWeight if present, else 1)
        double w = HAS_CENT_WEIGHT ? (double)centralityWeight : 1.0;

        bool haveReco = (reco_pt > (Float_t)RECO_DUMMY_CUT);
        if (!haveReco) continue;

        // reco quality cuts
        if (!(reco_area >= (Float_t)areaMin)) continue;
        if (!(reco_neutral_fraction <= (Float_t)CUT_NEUTRAL_FRACTION)) continue;

        // trigger match
        if (!(reco_trigger_match == kTRUE)) continue;

        // fill spectra per ptlead threshold
        for (int it = 0; it < N_LEAD; ++it) {
          if (reco_pt_lead >= (Float_t)PTLEAD_THR[it]) {
            hSpec[it]->Fill((double)reco_pt_corr, w);
          }
        }
      }

      // ---------------- write hists ----------------
      outC->cd();
      for (int it = 0; it < N_LEAD; ++it) hSpec[it]->Write();

      // ---------------- make PDF: overlay of all ptlead ----------------
      {
        TCanvas* c = new TCanvas("cSpec","cSpec",800,700);
        c->SetLogy();

        // base draw
        hSpec[0]->SetTitle("");
        hSpec[0]->GetXaxis()->SetTitle("#it{p}_{T, reco}^{corr} (GeV/#it{c})");
        hSpec[0]->GetYaxis()->SetTitle(HAS_CENT_WEIGHT ? "Weighted counts" : "Counts");
        hSpec[0]->GetXaxis()->SetRangeUser(-20, 60);
        hSpec[0]->Draw("E1");

        // styling (match your original)
        hSpec[0]->SetMarkerStyle(20);
        hSpec[0]->SetLineColor(kBlack);
        hSpec[0]->SetMarkerColor(kBlack);

        hSpec[1]->SetMarkerStyle(21);
        hSpec[1]->SetLineColor(kRed+1);
        hSpec[1]->SetMarkerColor(kRed+1);

        hSpec[2]->SetMarkerStyle(22);
        hSpec[2]->SetLineColor(kBlue+1);
        hSpec[2]->SetMarkerColor(kBlue+1);

        hSpec[3]->SetMarkerStyle(23);
        hSpec[3]->SetLineColor(kGreen+2);
        hSpec[3]->SetMarkerColor(kGreen+2);

        hSpec[1]->Draw("E1 SAME");
        hSpec[2]->Draw("E1 SAME");
        hSpec[3]->Draw("E1 SAME");

        TLegend* leg = new TLegend(0.60, 0.68, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.04);
        for (int it=0; it<N_LEAD; ++it) {
          leg->AddEntry(hSpec[it],
                        Form("#it{p}_{T}^{lead} #geq %.0f GeV", PTLEAD_THR[it]),
                        "lep");
        }
        leg->Draw();

        double Rval = 0.0;
        sscanf(rname.c_str(), "R%lf", &Rval);

        TString centLabel = cname.c_str();
        centLabel.ReplaceAll("CENT_", "");
        centLabel.ReplaceAll("MID_",  "");
        centLabel.ReplaceAll("PERI_", "");
        centLabel.ReplaceAll("_", "-");
        centLabel += " %";

        TString line1;
        line1.Form("R = %.1f, %s", Rval, centLabel.Data());

        TLatex tex;
        tex.SetNDC();
        tex.SetTextFont(42);
        tex.SetTextSize(0.035);

        tex.DrawLatex(0.60, 0.60, line1.Data());
        tex.DrawLatex(0.60, 0.555, "Au+Au  #sqrt{s_{NN}} = 200 GeV");
        tex.DrawLatex(0.60, 0.51, "THIS THESIS");

        TString pdfName;
        pdfName.Form("%s/spectra_ptlead_%s_%s.pdf", pdfCDir.Data(), rname.c_str(), cname.c_str());
        c->SaveAs(pdfName.Data());

        delete leg;
        delete c;
      }

      // cleanup
      for (int it = 0; it < N_LEAD; ++it) delete hSpec[it];
    }
  }

  fout->Write();
  fout->Close();
  fin->Close();

  std::cout << "Done. ROOT output: " << outfile << "\n";
  std::cout << "PDFs written under: " << pdfTop << "/R*/CENT_*/\n";
}
