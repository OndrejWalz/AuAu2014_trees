// compare_bayes_svd.C
//
// Compares unfolded spectra from Bayesian and SVD unfolding.
// Reads from two separate files produced by unfold_embedding.C:
//   out_embedding_Bayes/responses_embedding.root  (method="Bayes")
//   out_embedding_SVD/responses_embedding.root    (method="SVD")
//
// Expected directory structure in each file:
//   R0.2_CENT_0_10_ptlead0/
//     Unfolded_iter1, Unfolded_iter4, Unfolded_iter5, Unfolded_iter7  (Bayes file)
//     Unfolded_SVD_iter3, Unfolded_SVD_iter4, Unfolded_SVD_iter5      (SVD file)
//     hTrueTest, hMeasTest, ...
//
// Usage (ROOT prompt):
//   .x compare_bayes_svd.C(
//        "out_embedding_Bayes/responses_embedding.root",
//        "out_embedding_SVD/responses_embedding.root",
//        "comparison_output")

#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TDirectory.h"

#include <vector>
#include <string>
#include <iostream>

using std::string;
using std::vector;
using std::cout;

// ----------------------- config (must match unfold_embedding.C) ---------------

static const double kPtLeadCuts[] = {0.0, 5.0, 7.0, 9.0};
static const int    kNPtLeadCuts  = sizeof(kPtLeadCuts)/sizeof(kPtLeadCuts[0]);

static const vector<string> kCentralities = {"CENT_0_10", "MID_20_40", "PERI_60_80"};
static const vector<string> kRadii        = {"R0.2", "R0.3", "R0.4"};

// Bayes iterations — must match kBayesIters[] in unfold_embedding.C
static const int kBayesIters[] = {1, 4, 5, 7};
static const int kNBayesIters  = sizeof(kBayesIters)/sizeof(kBayesIters[0]);

// SVD regularization values — must match kRegValues[] in unfold_embedding.C
static const int kSVDRegs[]  = {3, 4, 5};
static const int kNSVDRegs   = sizeof(kSVDRegs)/sizeof(kSVDRegs[0]);

// Which Bayes iteration and SVD k to use for the bottom-pad direct ratio (0-based index).
static const int kBestBayesIdx = kNBayesIters - 1;  // iter 7
static const int kBestSVDIdx   = kNSVDRegs    - 1;  // k = 5

// Truth binning
static const int    nbins_truth = 10;
static const double bin_truth_edges[nbins_truth+1] = {
  0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 60
};

// ----------------------- helpers ----------------------------------------------

static void EnsureDir(const string& path) {
  if (gSystem->AccessPathName(path.c_str()))
    gSystem->mkdir(path.c_str(), true);
}

static TString NiceCentLabel(const std::string& tok) {
  int a = -1, b = -1;
  if (sscanf(tok.c_str(), "CENT_%d_%d", &a, &b) == 2 ||
      sscanf(tok.c_str(), "MID_%d_%d",  &a, &b) == 2 ||
      sscanf(tok.c_str(), "PERI_%d_%d", &a, &b) == 2)
    return Form("%d#font[52]{#minus}%d %%", a, b);
  return TString(tok.c_str());
}

static TH1D* GetHist(TFile* f, const string& dirName, const string& histName) {
  TDirectory* d = dynamic_cast<TDirectory*>(f->Get(dirName.c_str()));
  if (!d) {
    cout << "[warn] directory not found in " << f->GetName() << ": " << dirName << "\n";
    return nullptr;
  }
  TH1D* h = dynamic_cast<TH1D*>(d->Get(histName.c_str()));
  if (!h)
    cout << "[warn] histogram not found: " << dirName << "/" << histName
         << " in " << f->GetName() << "\n";
  return h;
}

// Clone + rebin to truth binning. Caller owns the result.
static TH1D* RebinToTruth(TH1D* hin, const string& newName) {
  if (!hin) return nullptr;
  if (hin->GetNbinsX() == nbins_truth) {
    TH1D* c = (TH1D*)hin->Clone(newName.c_str());
    c->SetDirectory(nullptr);
    return c;
  }
  TH1D* h = (TH1D*)hin->Rebin(nbins_truth, newName.c_str(), bin_truth_edges);
  h->SetDirectory(nullptr);
  return h;
}

static void DrawRefLines(double xlo, double xhi) {
  struct { double y; int col; int style; } L[] = {
    {1.00, kBlack,  1},
    {1.05, kGray,   2}, {0.95, kGray,   2},
    {1.10, kGray+2, 2}, {0.90, kGray+2, 2},
  };
  for (int i = 0; i < 5; ++i) {
    TLine* ln = new TLine(xlo, L[i].y, xhi, L[i].y);
    ln->SetLineColor(L[i].col);
    ln->SetLineStyle(L[i].style);
    ln->Draw();
  }
}

// ----------------------- main -------------------------------------------------

void compare_bayes_svd(const char* bayesFile,
                       const char* svdFile,
                       const char* outDir)
{
  gStyle->SetOptStat(0);
  EnsureDir(outDir);
  TH1::SetDefaultSumw2(kTRUE);

  TFile* fB = TFile::Open(bayesFile, "READ");
  if (!fB || fB->IsZombie()) {
    cout << "[error] Cannot open Bayes file: " << bayesFile << "\n"; return;
  }
  TFile* fS = TFile::Open(svdFile, "READ");
  if (!fS || fS->IsZombie()) {
    cout << "[error] Cannot open SVD file: " << svdFile << "\n";
    fB->Close(); delete fB; return;
  }

  cout << "Bayes file : " << bayesFile << "\n";
  cout << "SVD file   : " << svdFile   << "\n";
  cout << "Output dir : " << outDir    << "\n\n";

  // Palettes — Bayes: warm/filled, SVD: cool/open
  static const int kBayesCols[]  = {kRed-4, kOrange+7, kRed+2, kMagenta+1};
  static const int kBayesMarks[] = {20, 33, 22, 21};
  static const int kSVDCols[]    = {kAzure+2, kCyan+2, kGreen+2};
  static const int kSVDMarks[]   = {24, 26, 25};

  for (size_t iR = 0; iR < kRadii.size(); ++iR) {
    const string R = kRadii[iR];
    double Rval = 0.0;
    sscanf(R.c_str(), "R%lf", &Rval);

    for (size_t iC = 0; iC < kCentralities.size(); ++iC) {
      const string C = kCentralities[iC];

      for (int ic = 0; ic < kNPtLeadCuts; ++ic) {
        const double cut = kPtLeadCuts[ic];
        const string tag = R + "_" + C + Form("_ptlead%.0f", cut);
        cout << "=== " << tag << " ===\n";

        // ---- Bayes histos (from Bayes file) ----
        vector<TH1D*> hBayes(kNBayesIters, nullptr);
        for (int ib = 0; ib < kNBayesIters; ++ib) {
          TH1D* raw = GetHist(fB, tag, Form("Unfolded_iter%d", kBayesIters[ib]));
          hBayes[ib] = RebinToTruth(raw,
              Form("hBayes_iter%d_%s", kBayesIters[ib], tag.c_str()));
        }

        // ---- SVD histos (from SVD file) ----
        vector<TH1D*> hSVD(kNSVDRegs, nullptr);
        for (int is = 0; is < kNSVDRegs; ++is) {
          TH1D* raw = GetHist(fS, tag, Form("Unfolded_SVD_iter%d", kSVDRegs[is]));
          hSVD[is] = RebinToTruth(raw,
              Form("hSVD_reg%d_%s", kSVDRegs[is], tag.c_str()));
        }

        // ---- truth: prefer Bayes file, fall back to SVD file ----
        TH1D* hTruthRaw = GetHist(fB, tag, "hTrueTest");
        if (!hTruthRaw) hTruthRaw = GetHist(fS, tag, "hTrueTest");
        if (!hTruthRaw) {
          cout << "[warn] No hTrueTest for " << tag << ", skipping.\n";
          continue;
        }
        TH1D* hTruth = RebinToTruth(hTruthRaw, ("hTruth_"+tag).c_str());

        // ---------------------------------------------------------------
        // Canvas: 3 pads
        //   pSpec  (top,  45-100%) : all spectra, log-y
        //   pRatio (mid,  18-45%)  : every (Unfolded / Truth)
        //   pComp  (bot,   0-18%)  : Bayes(best) / SVD(best)
        // ---------------------------------------------------------------
        TCanvas* c = new TCanvas(("c_"+tag).c_str(), "", 800, 1100);

        TPad* pSpec  = new TPad(("pSpec_"+tag).c_str(),  "", 0.0, 0.45, 1.0, 1.00);
        TPad* pRatio = new TPad(("pRatio_"+tag).c_str(), "", 0.0, 0.18, 1.0, 0.45);
        TPad* pComp  = new TPad(("pComp_"+tag).c_str(),  "", 0.0, 0.00, 1.0, 0.18);

        pSpec ->SetLeftMargin(0.12); pSpec ->SetRightMargin(0.03);
        pSpec ->SetTopMargin (0.04); pSpec ->SetBottomMargin(0.01);
        pRatio->SetLeftMargin(0.12); pRatio->SetRightMargin(0.03);
        pRatio->SetTopMargin (0.01); pRatio->SetBottomMargin(0.01);
        pComp ->SetLeftMargin(0.12); pComp ->SetRightMargin(0.03);
        pComp ->SetTopMargin (0.01); pComp ->SetBottomMargin(0.35);

        pSpec->Draw(); pRatio->Draw(); pComp->Draw();

        // ==================== PAD 1: spectra ====================
        pSpec->cd();
        pSpec->SetLogy();

        hTruth->SetMarkerStyle(20);
        hTruth->SetMarkerColor(kBlack);
        hTruth->SetLineColor(kBlack);
        hTruth->SetMarkerSize(1.1);
        hTruth->GetXaxis()->SetLabelSize(0);
        hTruth->GetXaxis()->SetTitleSize(0);
        hTruth->GetYaxis()->SetTitle("dN/dp_{T}");
        hTruth->GetYaxis()->SetTitleOffset(1.25);
        hTruth->GetYaxis()->SetTitleSize(0.045);
        hTruth->GetYaxis()->SetLabelSize(0.040);
        hTruth->Draw("E1");

        TLegend* leg = new TLegend(0.55, 0.52, 0.94, 0.94);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
        leg->AddEntry(hTruth, "Truth (test)", "lp");

        for (int ib = 0; ib < kNBayesIters; ++ib) {
          if (!hBayes[ib]) continue;
          hBayes[ib]->SetMarkerStyle(kBayesMarks[ib]);
          hBayes[ib]->SetMarkerColor(kBayesCols[ib]);
          hBayes[ib]->SetLineColor(kBayesCols[ib]);
          hBayes[ib]->Draw("E1 SAME");
          leg->AddEntry(hBayes[ib], Form("Bayes %d it.", kBayesIters[ib]), "lp");
        }
        for (int is = 0; is < kNSVDRegs; ++is) {
          if (!hSVD[is]) continue;
          hSVD[is]->SetMarkerStyle(kSVDMarks[is]);
          hSVD[is]->SetMarkerColor(kSVDCols[is]);
          hSVD[is]->SetLineColor(kSVDCols[is]);
          hSVD[is]->Draw("E1 SAME");
          leg->AddEntry(hSVD[is], Form("SVD k=%d", kSVDRegs[is]), "lp");
        }
        leg->Draw();

        {
          TLatex lat;
          lat.SetNDC(true); lat.SetTextFont(42); lat.SetTextSize(0.038);
          lat.DrawLatex(0.16, 0.30, "Au+Au  #sqrt{#it{s}_{NN}} = 200 GeV");
          lat.DrawLatex(0.16, 0.23,
              Form("#it{R} = %.1f,  %s", Rval, NiceCentLabel(C).Data()));
          lat.DrawLatex(0.16, 0.16,
              Form("#it{p}_{T}^{lead} #geq %.0f GeV/#it{c}", cut));
        }

        // ==================== PAD 2: Unfolded / Truth ====================
        pRatio->cd();
        //pRatio->SetGridy();

        TH1D* firstR = nullptr;
        auto DrawRatio = [&](TH1D* hunf, int col, int mark) {
          if (!hunf) return;
          TH1D* r = (TH1D*)hunf->Clone(Form("r_%s", hunf->GetName()));
          r->SetDirectory(nullptr);
          r->Divide(hTruth);
          r->SetMarkerStyle(mark);
          r->SetMarkerColor(col);
          r->SetLineColor(col);
          if (!firstR) {
            r->SetTitle("");
            r->GetYaxis()->SetTitle("Unfolded / Truth");
            r->GetYaxis()->SetRangeUser(0.5, 1.5);
            r->GetYaxis()->SetTitleSize(0.075);
            r->GetYaxis()->SetLabelSize(0.065);
            r->GetYaxis()->SetTitleOffset(0.75);
            r->GetXaxis()->SetLabelSize(0);
            r->GetXaxis()->SetTitleSize(0);
            r->Draw("E1");
            firstR = r;
          } else {
            r->Draw("E1 SAME");
          }
        };

        for (int ib = 0; ib < kNBayesIters; ++ib)
          DrawRatio(hBayes[ib], kBayesCols[ib], kBayesMarks[ib]);
        for (int is = 0; is < kNSVDRegs; ++is)
          DrawRatio(hSVD[is], kSVDCols[is], kSVDMarks[is]);

        if (firstR)
          DrawRefLines(firstR->GetXaxis()->GetXmin(),
                       firstR->GetXaxis()->GetXmax());

        // ==================== PAD 3: Bayes(best) / SVD(best) ====================
        pComp->cd();
        //pComp->SetGridy();

        TH1D* hB = hBayes[kBestBayesIdx];
        TH1D* hS = hSVD  [kBestSVDIdx];

        if (hB && hS) {
          TH1D* hComp = (TH1D*)hB->Clone(("BayesSVD_"+tag).c_str());
          hComp->SetDirectory(nullptr);
          hComp->Divide(hS);

          hComp->SetMarkerStyle(20);
          hComp->SetMarkerColor(kViolet+2);
          hComp->SetLineColor(kViolet+2);
          hComp->SetTitle("");

          hComp->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
          hComp->GetXaxis()->SetTitleSize(0.13);
          hComp->GetXaxis()->SetLabelSize(0.11);
          hComp->GetXaxis()->SetTitleOffset(1.05);

          hComp->GetYaxis()->SetTitle(
              Form("Bayes %dit / SVD k=%d",
                   kBayesIters[kBestBayesIdx], kSVDRegs[kBestSVDIdx]));
          hComp->GetYaxis()->SetRangeUser(0.7, 1.3);
          hComp->GetYaxis()->SetTitleSize(0.095);
          hComp->GetYaxis()->SetLabelSize(0.090);
          hComp->GetYaxis()->SetTitleOffset(0.58);
          hComp->GetYaxis()->SetNdivisions(505);

          hComp->Draw("E1");
          DrawRefLines(hComp->GetXaxis()->GetXmin(),
                       hComp->GetXaxis()->GetXmax());
        } else {
          cout << "[warn] Missing best Bayes or SVD histo for " << tag
               << " — bottom pad empty.\n";
        }

        // ---------- save ----------
        const string base = string(outDir) + "/BayesVsSVD_" + tag;
        c->SaveAs((base + ".pdf").c_str());
        c->SaveAs((base + ".png").c_str());
        cout << "  saved " << base << ".png\n";

        // ---------- cleanup ----------
        delete c;
        delete leg;
        delete hTruth;
        for (int ib = 0; ib < kNBayesIters; ++ib) delete hBayes[ib];
        for (int is = 0; is < kNSVDRegs;    ++is) delete hSVD[is];

      } // ptlead
    } // centrality
  } // radius

  fB->Close(); delete fB;
  fS->Close(); delete fS;
  cout << "\nDone. Plots written to: " << outDir << "\n";
}