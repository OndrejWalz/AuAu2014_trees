// make_ratios.C
// ROOT5 / old C++ compatible (no nullptr, no C++11)
//
// Usage:
//   root -l -b -q 'make_ratios.C("hists.root","hists_withVeto.root","ratio_pdf","ratios.root")'

#include "TFile.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TKey.h"
#include "TClass.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"

#include <string>
#include <iostream>
#include <algorithm>

// ------------------------------------------------------------

bool Is1DHist(TObject* obj)
{
  if (!obj) return false;
  if (!obj->InheritsFrom(TH1::Class())) return false;
  if (obj->InheritsFrom(TH2::Class())) return false;

  TH1* h = dynamic_cast<TH1*>(obj);
  if (!h) return false;
  return (h->GetDimension() == 1);
}

void Nice1D(TH1* h, int mstyle, int lcolor)
{
  if (!h) return;
  h->SetMarkerStyle(mstyle);
  h->SetMarkerSize(0.9);
  h->SetLineColor(lcolor);
  h->SetMarkerColor(lcolor);
  h->SetLineWidth(2);
}

// ------------------------------------------------------------

void DrawRatioCanvas(TH1* hNomIn, TH1* hVarIn,
                     const char* outPdf,
                     const char* titleLine,
                     double ratioYmin,
                     double ratioYmax)
{
  if (!hNomIn || !hVarIn) return;

  TH1* hNom = (TH1*)hNomIn->Clone("hNom_clone");
  TH1* hVar = (TH1*)hVarIn->Clone("hVar_clone");
  hNom->SetDirectory(0);
  hVar->SetDirectory(0);

  if (hNom->GetNbinsX() != hVar->GetNbinsX() ||
      hNom->GetXaxis()->GetXmin() != hVar->GetXaxis()->GetXmin() ||
      hNom->GetXaxis()->GetXmax() != hVar->GetXaxis()->GetXmax()) {
    std::cerr << "  [SKIP] binning mismatch: " << outPdf << std::endl;
    delete hNom;
    delete hVar;
    return;
  }

  TH1* hRat = (TH1*)hVar->Clone("hRatio");
  hRat->SetDirectory(0);
  hRat->Divide(hNom);

  TCanvas* c = new TCanvas("cRatio","cRatio",900,850);

  TPad* p1 = new TPad("p1","p1",0.0,0.30,1.0,1.0);
  TPad* p2 = new TPad("p2","p2",0.0,0.00,1.0,0.30);
  p1->SetBottomMargin(0.02);
  p2->SetTopMargin(0.02);
  p2->SetBottomMargin(0.30);
  p1->Draw();
  p2->Draw();

  // --- top pad
  p1->cd();
  p1->SetTicks(1,1);

  double maxY = std::max(hNom->GetMaximum(), hVar->GetMaximum());
  if (maxY <= 0) maxY = 1.0;

  hNom->SetMaximum(1.35*maxY);
  hNom->SetMinimum(0.0);
  hNom->GetXaxis()->SetLabelSize(0);
  hNom->GetXaxis()->SetTitleSize(0);

  Nice1D(hNom, 20, kBlack);
  Nice1D(hVar, 24, kRed+1);

  hNom->Draw("E1");
  hVar->Draw("E1 SAME");

  TLegend* leg = new TLegend(0.60,0.73,0.88,0.88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.040);
  leg->AddEntry(hNom, "nominal", "lep");
  leg->AddEntry(hVar, "with veto", "lep");
  leg->Draw();

  TLatex tex;
  tex.SetNDC();
  tex.SetTextFont(42);
  tex.SetTextSize(0.042);
  tex.DrawLatex(0.14, 0.92, titleLine);

  // --- bottom pad
  p2->cd();
  p2->SetTicks(1,1);

  Nice1D(hRat, 20, kBlack);

  hRat->GetYaxis()->SetTitle("veto / nominal");
  hRat->GetYaxis()->CenterTitle();
  hRat->GetYaxis()->SetNdivisions(505);
  hRat->GetYaxis()->SetTitleSize(0.10);
  hRat->GetYaxis()->SetLabelSize(0.09);
  hRat->GetYaxis()->SetTitleOffset(0.55);

  hRat->GetXaxis()->SetTitle(hNom->GetXaxis()->GetTitle());
  hRat->GetXaxis()->SetTitleSize(0.12);
  hRat->GetXaxis()->SetLabelSize(0.10);

  hRat->SetMinimum(ratioYmin);
  hRat->SetMaximum(ratioYmax);
  hRat->Draw("E1");

  TLine* l1 = new TLine(hRat->GetXaxis()->GetXmin(),1.0,
                        hRat->GetXaxis()->GetXmax(),1.0);
  l1->SetLineStyle(2);
  l1->Draw();

  c->SaveAs(outPdf);

  delete l1;
  delete leg;
  delete c;
  delete hNom;
  delete hVar;
  delete hRat;
}

// ------------------------------------------------------------

void WalkDirsAndMakeRatios(TDirectory* dNom, TDirectory* dVar,
                           const std::string& outBase,
                           const std::string& relPath,
                           TFile* foutRatios)
{
  if (!dNom || !dVar) return;

  std::string outDir = outBase;
  if (!relPath.empty()) outDir += "/" + relPath;
  gSystem->mkdir(outDir.c_str(), kTRUE);

  TIter nextKey(dNom->GetListOfKeys());
  TKey* key = 0;

  while ((key = (TKey*)nextKey())) {
    const char* name = key->GetName();

    if (std::string(key->GetClassName()) == "TDirectoryFile") {
      TDirectory* subNom = (TDirectory*)dNom->Get(name);
      TDirectory* subVar = (TDirectory*)dVar->Get(name);
      if (!subNom || !subVar) continue;

      std::string nextRel = relPath.empty() ? name : relPath + "/" + name;
      WalkDirsAndMakeRatios(subNom, subVar, outBase, nextRel, foutRatios);
      continue;
    }

    TObject* oNom = dNom->Get(name);
    TObject* oVar = dVar->Get(name);
    if (!oNom || !oVar) continue;

    if (!Is1DHist(oNom) || !Is1DHist(oVar)) continue;

    TH1* hNom = (TH1*)oNom;
    TH1* hVar = (TH1*)oVar;

    std::string pdfName = outDir + "/" + name + std::string("_ratio.pdf");
    std::string title   = relPath.empty() ? name : relPath + "/" + name;

    std::cout << "Ratio: " << title << std::endl;

    DrawRatioCanvas(hNom, hVar, pdfName.c_str(), title.c_str(), 0.8, 1.2);

    if (foutRatios) {
      foutRatios->cd();
      TDirectory* outD = foutRatios;

      if (!relPath.empty()) {
        size_t pos = 0;
        while (pos < relPath.size()) {
          size_t next = relPath.find('/', pos);
          std::string part = relPath.substr(pos, next-pos);
          if (!outD->GetDirectory(part.c_str()))
            outD = outD->mkdir(part.c_str());
          else
            outD = outD->GetDirectory(part.c_str());
          if (next == std::string::npos) break;
          pos = next + 1;
        }
      }

      outD->cd();
      TH1* hRat = (TH1*)hVar->Clone(Form("%s_ratio", name));
      hRat->Divide(hNom);
      hRat->Write();
      delete hRat;
    }
  }
}

// ------------------------------------------------------------

void make_ratios(const char* nomFile,
                 const char* varFile,
                 const char* outPdfDir,
                 const char* outRoot)
{
  TH1::SetDefaultSumw2(kTRUE);
  gStyle->SetOptStat(0);
  gROOT->SetBatch(kTRUE);

  TFile* fNom = TFile::Open(nomFile,"READ");
  TFile* fVar = TFile::Open(varFile,"READ");
  if (!fNom || !fVar) return;

  gSystem->mkdir(outPdfDir, kTRUE);

  TFile* fOut = 0;
  if (outRoot && strlen(outRoot))
    fOut = TFile::Open(outRoot,"RECREATE");

  WalkDirsAndMakeRatios(fNom, fVar, outPdfDir, "", fOut);

  if (fOut) { fOut->Write(); fOut->Close(); }
  fNom->Close();
  fVar->Close();

  std::cout << "Done." << std::endl;
}
