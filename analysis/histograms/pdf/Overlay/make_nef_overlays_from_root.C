// make_nef_overlays_from_root.C
// ROOT5-safe: clone + SetDirectory(0) before drawing.
// Output PDFs are written to the CURRENT directory (no new folders created).
//
// Updates vs previous version:
//  - Increased axis label/title sizes
//  - Increased pad margins to avoid clipping
//  - Added vertical dashed red line at x = 0.95

#include "TFile.h"
#include "TDirectoryFile.h"
#include "TKey.h"
#include "TH1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TString.h"
#include "TLine.h"

#include <iostream>
#include <vector>
#include <string>

static void NormalizeProb(TH1* h)
{
  if (!h) return;
  double I = h->Integral();
  if (I > 0) h->Scale(1.0 / I);
}

static TString CentralityLabelFromDir(const std::string& cname)
{
  if (cname == "CENT_0_10")  return "0-10 %";
  if (cname == "MID_20_40")  return "20-40 %";
  if (cname == "PERI_60_80") return "60-80 %";
  return "";
}

static void DrawHeaderTopLeft(const TString& centLabel = "")
{
  if (!gPad) return;

  const double x    = gPad->GetLeftMargin() + 0.02;
  const double yTop = 1.0 - gPad->GetTopMargin() - 0.02;
  const double dy   = 0.050;

  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(13); // left-top

  tx.SetTextSize(0.040);
  if (centLabel.Length() > 0) {
    tx.DrawLatex(x, yTop,
                 Form("STAR Au+Au %s  #sqrt{s_{NN}} = 200 GeV", centLabel.Data()));
  } else {
    tx.DrawLatex(x, yTop, "STAR Au+Au  #sqrt{s_{NN}} = 200 GeV");
  }

  tx.SetTextSize(0.034);
  tx.DrawLatex(x, yTop - dy, "THIS THESIS");
}

// ---------- axis cosmetics ----------
static void ApplyAxisStyle1D(TH1* h)
{
  if (!h) return;

  const double lab = 0.045;
  const double tit = 0.050;

  if (h->GetXaxis()) {
    h->GetXaxis()->SetLabelSize(lab);
    h->GetXaxis()->SetTitleSize(tit);
    h->GetXaxis()->SetTitleOffset(1.10);
  }
  if (h->GetYaxis()) {
    h->GetYaxis()->SetLabelSize(lab);
    h->GetYaxis()->SetTitleSize(tit);
    h->GetYaxis()->SetTitleOffset(1.35);
  }
}

// ---------- vertical line at fixed x ----------
static void DrawVLineX(double x)
{
  if (!gPad) return;

  gPad->Update();

  // Use the drawn frame range (works after at least one Draw + Update)
  const double xmin = gPad->GetUxmin();
  const double xmax = gPad->GetUxmax();
  if (xmax <= xmin) return;

  // Clamp x to visible range (optional but safe)
  double xx = x;
  if (xx < xmin) xx = xmin;
  if (xx > xmax) xx = xmax;

  const double L = gPad->GetLeftMargin();
  const double R = gPad->GetRightMargin();
  const double B = gPad->GetBottomMargin();
  const double T = gPad->GetTopMargin();

  const double frac = (xx - xmin) / (xmax - xmin);
  const double xNDC = L + (1.0 - L - R) * frac;

  const double y1NDC = B;
  const double y2NDC = 1.0 - T;

  TLine* ln = new TLine(xNDC, y1NDC, xNDC, y2NDC);
  ln->SetNDC(kTRUE);
  ln->SetLineColor(kRed);
  ln->SetLineStyle(2);
  ln->SetLineWidth(3);
  ln->Draw("SAME");

  gPad->Modified();
  gPad->Update();
}


static double SmallestPositiveBin(TH1* h)
{
  double minPos = 0.0;
  if (!h) return minPos;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double c = h->GetBinContent(i);
    if (c > 0 && (minPos <= 0 || c < minPos)) minPos = c;
  }
  return minPos;
}

// Reads histogram "hNEF" from directory R*/CENT* in the input file.
// Returns a detached clone (caller must delete).
static TH1D* GetNEFClone(TFile* fin,
                         const char* rname,
                         const char* cname)
{
  if (!fin) return 0;

  TDirectoryFile* dirR = (TDirectoryFile*)fin->Get(rname);
  if (!dirR) return 0;

  TDirectoryFile* dirC = (TDirectoryFile*)dirR->Get(cname);
  if (!dirC) return 0;

  TH1* hInBase = (TH1*)dirC->Get("hNEF");
  if (!hInBase) {
    std::cout << "Missing hist: " << rname << "/" << cname << "/hNEF\n";
    return 0;
  }

  TH1D* hIn = dynamic_cast<TH1D*>(hInBase);
  if (!hIn) {
    std::cout << "Found hNEF but it is not TH1D, it is: "
              << hInBase->ClassName() << " in " << rname << "/" << cname << "\n";
    return 0;
  }

  TString newName = Form("hNEF_%s_%s_clone", rname, cname);
  TH1D* h = (TH1D*)hIn->Clone(newName.Data());
  if (!h) return 0;

  h->SetDirectory(0);
  h->SetTitle("");
  return h;
}

void make_nef_overlays_from_root(const char* infile = "nef_inputs.root")
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetTitle(0);

  // Prevent ROOT from auto-attaching histograms to directories
  TH1::AddDirectory(kFALSE);

  // Open input
  TFile* fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Error: cannot open input file: " << infile << "\n";
    std::cerr << "Tip: run ROOT from the directory where the file is, or pass a full path:\n";
    std::cerr << "  root -l -b -q 'make_nef_overlays_from_root.C(\"/full/path/nef_hists.root\")'\n";
    return;
  }

  // Determine centralities from R0.2 keys
  TDirectoryFile* dirR02 = (TDirectoryFile*)fin->Get("R0.2");
  if (!dirR02) {
    std::cerr << "Error: missing directory R0.2 in " << infile << "\n";
    fin->Close();
    delete fin;
    return;
  }

  std::vector<std::string> cents;
  TIter it(dirR02->GetListOfKeys());
  TKey* k = 0;
  while ((k = (TKey*)it())) {
    if (!k->IsFolder()) continue;
    cents.push_back(k->GetName());
  }

  std::cout << "Found " << cents.size() << " centralities:\n";
  for (size_t i = 0; i < cents.size(); ++i) std::cout << "  " << cents[i] << "\n";

  const char* Rnames[3] = {"R0.2", "R0.3", "R0.4"};
  const int   cols[3]   = {kBlack, kBlue+1, kGreen+2};

  // Make overlays in CURRENT directory
  for (size_t ic = 0; ic < cents.size(); ++ic) {
    const std::string& cname = cents[ic];
    const TString centLabel  = CentralityLabelFromDir(cname);

    TH1D* h[3] = {0,0,0};
    for (int ir = 0; ir < 3; ++ir) {
      h[ir] = GetNEFClone(fin, Rnames[ir], cname.c_str());
      if (h[ir]) NormalizeProb(h[ir]);
    }

    if (!h[0] && !h[1] && !h[2]) {
      for (int ir = 0; ir < 3; ++ir) delete h[ir];
      continue;
    }

    TCanvas* c = new TCanvas(Form("c_nef_%s", cname.c_str()), "", 800, 600);
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);
    c->SetTopMargin(0.06);
    c->SetRightMargin(0.05);
    c->SetLogy(1);
    c->SetTickx(0);
    c->SetTicky(0);

    TH1D* href = h[0] ? h[0] : (h[1] ? h[1] : h[2]);
    href->SetLineColor(cols[0]);
    href->SetLineWidth(2);
    href->GetXaxis()->SetTitle("NEF (-)");
    href->GetYaxis()->SetTitle("Probability (a.u.)");
    ApplyAxisStyle1D(href);

    double minPos = SmallestPositiveBin(href);
    if (minPos > 0) href->SetMinimum(0.5 * minPos);
    href->Draw("HIST");

    for (int ir = 0; ir < 3; ++ir) {
      if (!h[ir] || h[ir] == href) continue;
      h[ir]->SetLineColor(cols[ir]);
      h[ir]->SetLineWidth(2);
      h[ir]->Draw("HIST SAME");
    }

    // vertical cut line
    DrawVLineX(0.95);

    TLegend* leg = new TLegend(0.62, 0.22, 0.88, 0.40);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    for (int ir = 0; ir < 3; ++ir) {
      if (h[ir]) leg->AddEntry(h[ir], Rnames[ir], "l");
    }
    leg->Draw();

    DrawHeaderTopLeft(centLabel);

    TString outPdf = Form("NEF_overlay_%s.pdf", cname.c_str());
    std::cout << "Saving: " << outPdf.Data() << "\n";
    c->SaveAs(outPdf.Data());

    delete leg;
    delete c;
    for (int ir = 0; ir < 3; ++ir) delete h[ir];
  }

  fin->Close();
  delete fin;
}
