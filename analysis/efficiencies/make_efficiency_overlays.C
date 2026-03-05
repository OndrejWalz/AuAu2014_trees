// make_efficiency_overlays.C  (ROOT5/CINT-friendly)
// Usage:
//   root -l -b -q 'make_efficiency_overlays.C("efficiencies.root","pdf/EffOverlays")'

#include "TFile.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TROOT.h"
#include "TCollection.h"

#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cstring>

static void ApplyNiceStyle()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.12);

  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTextFont(42);

  gStyle->SetTitleSize(0.045, "XYZ");
  gStyle->SetLabelSize(0.040, "XYZ");

  gStyle->SetEndErrorSize(4);
  gStyle->SetLineWidth(2);
}

static bool HasPrefix(const std::string& s, const char* prefix)
{
  if (!prefix) return false;
  size_t lp = strlen(prefix);
  if (s.size() < lp) return false;
  for (size_t i = 0; i < lp; ++i) if (s[i] != prefix[i]) return false;
  return true;
}

static bool Contains(const std::string& s, const char* sub)
{
  if (!sub) return false;
  return (s.find(sub) != std::string::npos);
}

static void PushUniqueStr(std::vector<std::string>& v, const std::string& s)
{
  for (size_t i = 0; i < v.size(); ++i) if (v[i] == s) return;
  v.push_back(s);
}

static void PushUniqueInt(std::vector<int>& v, int x)
{
  for (size_t i = 0; i < v.size(); ++i) if (v[i] == x) return;
  v.push_back(x);
}

static bool ExtractRadiusToken(const std::string& tag, std::string& rOut)
{
  size_t us = tag.find('_');
  if (us == std::string::npos) return false;
  rOut = tag.substr(0, us);
  if (rOut.size() < 2) return false;
  if (rOut[0] != 'R') return false;
  return true;
}

// For any tag like "R0.2_CENT_0_10_ptlead5" or "R0.2_MID_20_40_ptlead5":
// extract centrality token between first "_" and "_ptlead"
static bool ExtractAnyCentToken(const std::string& tag, std::string& centOut)
{
  size_t us = tag.find('_');           // after R0.x
  size_t ip = tag.rfind("_ptlead");
  if (us == std::string::npos || ip == std::string::npos) return false;
  if (ip <= us + 1) return false;
  centOut = tag.substr(us + 1, ip - (us + 1)); // "CENT_0_10" or "MID_20_40" etc.
  return true;
}

static bool ExtractPtLead(const std::string& tag, int& ptleadOut)
{
  size_t ip = tag.rfind("_ptlead");
  if (ip == std::string::npos) return false;
  std::string after = tag.substr(ip + strlen("_ptlead"));
  if (after.size() == 0) return false;
  ptleadOut = atoi(after.c_str());
  return true;
}

static void SortRadii(std::vector<std::string>& radii)
{
  if (radii.size() < 2) return;
  for (size_t i = 0; i < radii.size(); ++i) {
    for (size_t j = i + 1; j < radii.size(); ++j) {
      double ri = atof(radii[i].c_str() + 1);
      double rj = atof(radii[j].c_str() + 1);
      if (rj < ri) {
        std::string tmp = radii[i];
        radii[i] = radii[j];
        radii[j] = tmp;
      }
    }
  }
}

static void SortStrings(std::vector<std::string>& v)
{
  if (v.size() < 2) return;
  for (size_t i = 0; i < v.size(); ++i) {
    for (size_t j = i + 1; j < v.size(); ++j) {
      if (v[j] < v[i]) {
        std::string tmp = v[i];
        v[i] = v[j];
        v[j] = tmp;
      }
    }
  }
}

static void SortInts(std::vector<int>& v)
{
  if (v.size() < 2) return;
  for (size_t i = 0; i < v.size(); ++i) {
    for (size_t j = i + 1; j < v.size(); ++j) {
      if (v[j] < v[i]) {
        int tmp = v[i];
        v[i] = v[j];
        v[j] = tmp;
      }
    }
  }
}

static TString PrettyCentrality(const std::string& centToken)
{
  int a = -1, b = -1;

  if (sscanf(centToken.c_str(), "CENT_%d_%d", &a, &b) == 2)
    return Form("%d#font[52]{#minus}%d %%", a, b);

  if (sscanf(centToken.c_str(), "MID_%d_%d", &a, &b) == 2)
    return Form("%d#font[52]{#minus}%d %%", a, b);

  if (sscanf(centToken.c_str(), "PERI_%d_%d", &a, &b) == 2)
    return Form("%d#font[52]{#minus}%d %%", a, b);

  return TString(centToken.c_str());
}

static void DrawTextBlock(double x, double yTop,
                          int align, double size,
                          const char* l1, const char* l2)
{
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextAlign(align);
  lat.SetTextSize(size);
  if (l1 && strlen(l1)) lat.DrawLatex(x, yTop, l1);
  if (l2 && strlen(l2)) lat.DrawLatex(x, yTop - 0.06, l2);
}

// Efficiency builder from num_reb/den_reb (fallback to fine eff)
static TH1D* GetEfficiencyHist(TDirectory* d, const TString& tag, const char* mode)
{
  if (!d) return 0;

  TString hDenReb, hNumReb;
  if (strcmp(mode, "match") == 0) {
    hDenReb.Form("h_match_mc_den_reb_%s", tag.Data());
    hNumReb.Form("h_match_mc_num_reb_%s", tag.Data());
  } else {
    hDenReb.Form("h_trig_den_reb_%s", tag.Data());
    hNumReb.Form("h_trig_num_reb_%s", tag.Data());
  }

  TH1D* hDen = (TH1D*)d->Get(hDenReb);
  TH1D* hNum = (TH1D*)d->Get(hNumReb);
  if (hDen && hNum) {
    TH1D* hEff = (TH1D*)hNum->Clone(Form("tmp_%s_eff_%s", mode, tag.Data()));
    hEff->SetDirectory(0);
    hEff->SetTitle("");
    hEff->Divide(hNum, hDen, 1.0, 1.0, "b");
    return hEff;
  }

  TString hEffFine = (strcmp(mode, "match") == 0) ? "h_match_mc_eff" : "h_trig_eff";
  TH1D* hf = (TH1D*)d->Get(hEffFine);
  if (hf) {
    TH1D* hc2 = (TH1D*)hf->Clone(Form("tmp_%s_efffine_%s", mode, tag.Data()));
    hc2->SetDirectory(0);
    return hc2;
  }

  return 0;
}

void make_efficiency_overlays(const char* infile = "efficiencies.root",
                              const char* outdir = "pdf/EffOverlays")
{
  ApplyNiceStyle();
  TH1::SetDefaultSumw2(kTRUE);
  gSystem->mkdir(outdir, kTRUE);

  TFile* f = TFile::Open(infile, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Error: cannot open " << infile << std::endl;
    return;
  }

  std::vector<std::string> radii;
  std::vector<std::string> cents;
  std::vector<int> ptleads;

  // ---- Discover all tags, including MID_/PERI_ ----
  TIter nextKey(f->GetListOfKeys());
  TKey* key = 0;
  while ((key = (TKey*)nextKey())) {
    if (std::string(key->GetClassName()) != "TDirectoryFile") continue;

    std::string dname = key->GetName();
    if (!HasPrefix(dname, "R")) continue;
    if (!Contains(dname, "_ptlead")) continue;

    // require some recognizable centrality token
    if (!Contains(dname, "CENT_") && !Contains(dname, "MID_") && !Contains(dname, "PERI_"))
      continue;

    std::string rTok, cTok;
    int pTok = -999;

    if (ExtractRadiusToken(dname, rTok)) PushUniqueStr(radii, rTok);
    if (ExtractAnyCentToken(dname, cTok)) PushUniqueStr(cents, cTok);
    if (ExtractPtLead(dname, pTok)) PushUniqueInt(ptleads, pTok);
  }

  if (radii.size() == 0 || cents.size() == 0 || ptleads.size() == 0) {
    std::cerr << "Error: no suitable directories found in " << infile << std::endl;
    f->Close();
    return;
  }

  SortRadii(radii);
  SortStrings(cents);
  SortInts(ptleads);

  const int colors[]  = {kBlue+1, kRed+1, kGreen+2, kMagenta+1, kOrange+7};
  const int markers[] = {20, 21, 22, 33, 29};
  const int nStyles = 5;

  const char* hdr1 = "Au+Au #sqrt{#it{s}_{NN}} = 200 GeV";
  const char* hdr2 = "THIS THESIS";

  for (size_t ic = 0; ic < cents.size(); ++ic) {
    std::string centTok = cents[ic];
    TString centPretty = PrettyCentrality(centTok);

    for (size_t ip = 0; ip < ptleads.size(); ++ip) {
      int ptlead = ptleads[ip];

      TString centLine; centLine.Form("%s", centPretty.Data());
      TString leadLine; leadLine.Form("#it{p}_{T}^{lead} #geq %d GeV/#it{c}", ptlead);

      // ===================== MATCHING (UPDATED LAYOUT) =====================
      {
        TCanvas* c = new TCanvas("c_match_overlay", "c_match_overlay", 800, 600);

        // Legend in bottom-right but shifted left to leave space for text beside it
        // (Also lifted a bit to avoid the x-axis label.)
        TLegend* leg = new TLegend(0.58, 0.16, 0.74, 0.34);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.040);

        bool drew = false;

        for (size_t ir2 = 0; ir2 < radii.size(); ++ir2) {
          std::string r = radii[ir2];
          TString tag; tag.Form("%s_%s_ptlead%d", r.c_str(), centTok.c_str(), ptlead);

          TDirectory* d = (TDirectory*)f->Get(tag);
          if (!d) continue;

          TH1D* hc = GetEfficiencyHist(d, tag, "match");
          if (!hc) continue;

          int st = (int)(ir2 % nStyles);
          hc->SetLineColor(colors[st]);
          hc->SetMarkerColor(colors[st]);
          hc->SetMarkerStyle(markers[st]);
          hc->SetMarkerSize(1.0);
          hc->SetLineWidth(2);
          hc->SetStats(0);

          hc->GetXaxis()->SetTitle("#it{p}_{T}^{MC} (GeV/#it{c})");
          hc->GetYaxis()->SetTitle("Matching efficiency");
          hc->GetYaxis()->SetRangeUser(0.0, 1.10);

          if (!drew) { hc->Draw("E1"); drew = true; }
          else       { hc->Draw("E1 SAME"); }

          leg->AddEntry(hc, r.c_str(), "lep");
        }

        if (drew) {
          // Header stays upper-left
          DrawTextBlock(0.16, 0.88, 13, 0.045, hdr1, hdr2);

          // Put centrality + ptlead next to legend (bottom-right group)
          // Align left-top so it sits to the right of the legend box.
          DrawTextBlock(0.76, 0.32, 13, 0.040, centLine.Data(), leadLine.Data());

          leg->Draw();

          TString out;
          out.Form("%s/MatchEff_overlay_%s_ptlead%d.pdf", outdir, centTok.c_str(), ptlead);
          c->SaveAs(out);
        }

        delete leg;
        delete c;
      }

      // ===================== TRIGGER (UNCHANGED) =====================
      {
        TCanvas* c = new TCanvas("c_trig_overlay", "c_trig_overlay", 800, 600);

        // Legend left, below the centrality/pTlead block
        TLegend* leg = new TLegend(0.18, 0.14, 0.42, 0.32);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.040);

        bool drew = false;
        double yMax = 0.0;

        // yMax scan
        for (size_t ir = 0; ir < radii.size(); ++ir) {
          std::string r = radii[ir];
          TString tag; tag.Form("%s_%s_ptlead%d", r.c_str(), centTok.c_str(), ptlead);
          TDirectory* d = (TDirectory*)f->Get(tag);
          if (!d) continue;
          TH1D* hEff = GetEfficiencyHist(d, tag, "trig");
          if (!hEff) continue;
          double m = hEff->GetMaximum();
          if (m > yMax) yMax = m;
          delete hEff;
        }
        if (yMax <= 0.0) yMax = 1.0;

        for (size_t ir2 = 0; ir2 < radii.size(); ++ir2) {
          std::string r = radii[ir2];
          TString tag; tag.Form("%s_%s_ptlead%d", r.c_str(), centTok.c_str(), ptlead);

          TDirectory* d = (TDirectory*)f->Get(tag);
          if (!d) continue;

          TH1D* hc = GetEfficiencyHist(d, tag, "trig");
          if (!hc) continue;

          int st = (int)(ir2 % nStyles);
          hc->SetLineColor(colors[st]);
          hc->SetMarkerColor(colors[st]);
          hc->SetMarkerStyle(markers[st]);
          hc->SetMarkerSize(1.0);
          hc->SetLineWidth(2);
          hc->SetStats(0);

          hc->GetXaxis()->SetTitle("#it{p}_{T}^{reco,corr} (GeV/#it{c})");
          hc->GetYaxis()->SetTitle("Trigger efficiency");
          hc->GetXaxis()->SetRangeUser(-20.0, 60.0);
          hc->GetYaxis()->SetRangeUser(0.0, 1.15 * yMax);

          if (!drew) { hc->Draw("E1"); drew = true; }
          else       { hc->Draw("E1 SAME"); }

          leg->AddEntry(hc, r.c_str(), "lep");
        }

        if (drew) {
          // Header upper-left
          DrawTextBlock(0.16, 0.88, 13, 0.045, hdr1, hdr2);

          // Centrality + ptlead: LEFT-MIDDLE under header
          DrawTextBlock(0.16, 0.74, 13, 0.040, centLine.Data(), leadLine.Data());

          leg->Draw();

          TString out;
          out.Form("%s/TrigEff_overlay_%s_ptlead%d.pdf", outdir, centTok.c_str(), ptlead);
          c->SaveAs(out);
        }

        delete leg;
        delete c;
      }
    }
  }

  f->Close();
  std::cout << "Done. PDFs written to: " << outdir << std::endl;
}
