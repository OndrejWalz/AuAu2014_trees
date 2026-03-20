// rcp_check.cxx
// Usage (inside ROOT):
//   .x rcp_check.cxx+(".../unfolded_data.root","out_rcp","hUnfoldedTruthBins_matchCorr")
//
// It loops radii + ptlead and computes R_CP = (CENT/NcollCENT)/(PERI/NcollPERI).

#include "TFile.h"
#include "TDirectory.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TStyle.h"

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

using std::string;
using std::vector;
using std::cout;
using std::endl;

// ------------------------- config -----------------------------

// Ncoll per centrality bin: 0-5,5-10,10-20,20-30,30-40,40-50,50-60,60-70,70-80
static const double Ncoll[9] = {
  1066.5, 852.8, 606.9, 375.9, 222.6, 124.0, 64.0, 30.6, 13.7
};

static const double NcollErr[9] = {
  27.8, 23.4, 30.6, 33.4, 30.3, 24.6, 17.8, 11.4, 6.2
};

// TEMPORARY: equal-weight average to approximate 0-10 and 60-80
static const double Ncoll_0_10  = 0.5*(Ncoll[0] + Ncoll[1]);
static const double NcollErr_0_10 =
  0.5*std::sqrt(NcollErr[0]*NcollErr[0] + NcollErr[1]*NcollErr[1]);

static const double Ncoll_60_80 = 0.5*(Ncoll[7] + Ncoll[8]);
static const double NcollErr_60_80 =
  0.5*std::sqrt(NcollErr[7]*NcollErr[7] + NcollErr[8]*NcollErr[8]);

static const vector<string> kRadii = {"R0.2", "R0.3", "R0.4"};
static const double kPtLeadCuts[] = {0.0, 5.0, 7.0, 9.0};
static const int    kNPtLeadCuts  = sizeof(kPtLeadCuts)/sizeof(kPtLeadCuts[0]);

static const string kCent  = "CENT_0_10";
static const string kPeri  = "PERI_60_80";

// ------------------------- helpers ----------------------------

static void EnsureDir(const string& path){
  if (gSystem->AccessPathName(path.c_str()))
    gSystem->mkdir(path.c_str(), /*recursive=*/true);
}

static TH1D* GetHistClone(TFile* f, const string& dirName, const char* hname, const char* cloneName)
{
  if (!f) return 0;
  TDirectory* d = (TDirectory*) f->Get(dirName.c_str());
  if (!d) return 0;

  TH1D* h = 0;
  d->GetObject(hname, h);
  if (!h) return 0;

  TH1D* hc = (TH1D*) h->Clone(cloneName);
  hc->SetDirectory(0);
  return hc;
}

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

static TH1D* MakeRcp(const TH1D* hCent, const TH1D* hPeri,
                     double ncollCent, double ncollPeri,
                     const char* name)
{
  if (!hCent || !hPeri) return 0;
  if (!SameBinning1D(hCent, hPeri, 1e-6)) {
    cout << "  [ERROR] binning mismatch between cent/peri for " << name << endl;
    return 0;
  }

  TH1D* rcp = (TH1D*) hCent->Clone(name);
  rcp->SetDirectory(0);
  rcp->Reset();

  // R_CP = (C/NcC)/(P/NcP) = (C/P)*(NcP/NcC)
  const double scale = (ncollCent > 0.0) ? (ncollPeri / ncollCent) : 0.0;

  for (int ib = 1; ib <= rcp->GetNbinsX(); ++ib) {
    const double C  = hCent->GetBinContent(ib);
    const double eC = hCent->GetBinError(ib);
    const double P  = hPeri->GetBinContent(ib);
    const double eP = hPeri->GetBinError(ib);

    if (P > 0.0 && ncollCent > 0.0 && ncollPeri > 0.0) {
      const double val = (C / P) * scale;

      double err = 0.0;
      if (C > 0.0) {
        const double rel2 = (eC*eC)/(C*C) + (eP*eP)/(P*P);
        err = val * std::sqrt(rel2);
      } else {
        // if C==0, propagate absolute error from eC only (still reasonable)
        err = (eC / P) * scale;
      }

      rcp->SetBinContent(ib, val);
      rcp->SetBinError(ib, err);
    } else {
      rcp->SetBinContent(ib, 0.0);
      rcp->SetBinError(ib, 0.0);
    }
  }

  return rcp;
}

static void SaveQuickPdf(const TH1D* h, const char* pdf,
                         double ncollC, double ncollP,
                         double Rval, double ptlead)
{
  if (!h) return;

  TCanvas* c = new TCanvas(Form("c_%s", h->GetName()), "c", 800, 650);
  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.12);
  c->SetTopMargin(0.05);
  c->SetRightMargin(0.04);

  TH1D* hc = (TH1D*) h->Clone("hc_tmp");
  hc->SetDirectory(0);
  hc->SetTitle("");
  hc->GetXaxis()->SetTitle("#it{p}_{T, jet} (GeV)");
  hc->GetYaxis()->SetTitle("R_{CP}");
  hc->GetYaxis()->SetRangeUser(-0.5, 2.0);

  hc->Draw("E1");

  TLine* line = new TLine(hc->GetXaxis()->GetXmin(), 1.0,
                          hc->GetXaxis()->GetXmax(), 1.0);
  line->SetLineStyle(2);
  line->Draw("same");

  TLatex* lat = new TLatex();
  lat->SetNDC();
  lat->SetTextSize(0.035);
  lat->DrawLatex(0.15, 0.85, Form("N_{coll}^{cent} = %.3f", ncollC));
  lat->DrawLatex(0.15, 0.80, Form("N_{coll}^{peri} = %.3f", ncollP));

  lat->SetTextAlign(33); // right-top
  lat->DrawLatex(0.92, 0.90, "Au+Au  #sqrt{#it{s}_{NN}} = 200 GeV");
  lat->DrawLatex(0.92, 0.84, Form("#it{p}_{T}^{lead} #geq %.0f GeV/#it{c}", ptlead));
  lat->DrawLatex(0.92, 0.78, Form("#it{R} = %.1f", Rval));

  c->SaveAs(pdf);
  {
    TString png = pdf;
    png.ReplaceAll(".pdf", ".png");
    c->SaveAs(png);
  }

  delete lat;
  delete line;
  delete hc;
  delete c;
}

// ------------------------- main ------------------------------

void rcp_check(const char* unfoldedFile,
               const char* outDir = "out_rcp",
               const char* spectrumHistName = "hUnfoldedTruthBins_matchCorr")
{
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2(kTRUE);

  EnsureDir(outDir);
  EnsureDir(string(outDir) + "/pdf");

  TFile* fin = TFile::Open(unfoldedFile, "READ");
  if (!fin || fin->IsZombie()) {
    cout << "[error] cannot open input: " << unfoldedFile << endl;
    return;
  }

  const string outRoot = string(outDir) + "/rcp.root";
  TFile* fout = TFile::Open(outRoot.c_str(), "RECREATE");
  if (!fout || fout->IsZombie()) {
    cout << "[error] cannot create output: " << outRoot << endl;
    fin->Close(); delete fin;
    return;
  }

  cout << "\n=== R_CP quick check ===\n";
  cout << "Input : " << unfoldedFile << "\n";
  cout << "Hist  : " << spectrumHistName << "\n";
  cout << "Out   : " << outRoot << "\n\n";

  cout << "Using Ncoll(0-10)  = " << Ncoll_0_10  << " +- " << NcollErr_0_10  << "\n";
  cout << "Using Ncoll(60-80) = " << Ncoll_60_80 << " +- " << NcollErr_60_80 << "\n";
  cout << "Scale factor (peri/cent) = " << (Ncoll_60_80 / Ncoll_0_10) << "\n\n";

  // Loop radii + ptlead
  for (size_t iR = 0; iR < kRadii.size(); ++iR) {
    const string R = kRadii[iR];

    for (int ip = 0; ip < kNPtLeadCuts; ++ip) {
      const double cut = kPtLeadCuts[ip];

      const string tagCent = R + "_" + kCent + Form("_ptlead%.0f", cut);
      const string tagPeri = R + "_" + kPeri + Form("_ptlead%.0f", cut);

      TH1D* hC = GetHistClone(fin, tagCent, spectrumHistName,
                              Form("hCent_%s_ptlead%.0f", R.c_str(), cut));
      TH1D* hP = GetHistClone(fin, tagPeri, spectrumHistName,
                              Form("hPeri_%s_ptlead%.0f", R.c_str(), cut));

      if (!hC || !hP) {
        cout << "  [warn] missing hist for " << R << " ptlead " << cut
             << " (cent ok? " << (hC!=0) << ", peri ok? " << (hP!=0) << ")\n";
        if (hC) delete hC;
        if (hP) delete hP;
        continue;
      }

      const double nC = Ncoll_0_10;
      const double nP = Ncoll_60_80;

      TH1D* rcp = MakeRcp(hC, hP, nC, nP,
                          Form("hRcp_%s_ptlead%.0f", R.c_str(), cut));
      if (!rcp) {
        delete hC; delete hP;
        continue;
      }

      // write into directory per tag (optional structure)
      const string outTag = R + Form("_ptlead%.0f", cut);
      TDirectory* d = fout->mkdir(outTag.c_str());
      if (!d) d = fout->GetDirectory(outTag.c_str());
      d->cd();

      hC->Write("hCent", TObject::kOverwrite);
      hP->Write("hPeri", TObject::kOverwrite);
      rcp->Write("hRcp", TObject::kOverwrite);

      // quick pdf
      TString pdf = Form("%s/pdf/RCP_%s_ptlead%.0f.pdf", outDir, R.c_str(), cut);
      double Rval = 0.0;
      if (sscanf(R.c_str(), "R%lf", &Rval) != 1) Rval = 0.0;
      SaveQuickPdf(rcp, pdf.Data(), nC, nP, Rval, cut);

      delete hC;
      delete hP;
      delete rcp;
    }
  }

  fout->Write();
  fout->Close();
  delete fout;

  fin->Close();
  delete fin;

  cout << "\nDone. Wrote: " << outRoot << "\n";
}
