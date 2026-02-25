// make_qa_pdfs.C
//
// Produces QA PDFs + Jet PDFs.
// Adds a vertical red dashed line at the jet-area cut in:
//   (a) 1D area histograms
//   (b) 2D reco_pt_corr vs reco_area histograms
//
// NOTE: This does NOT apply the cut. It only draws the line.

#include "TFile.h"
#include "TDirectoryFile.h"
#include "TKey.h"
#include "TTree.h"
#include "TObject.h"
#include "TList.h"

#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TString.h"
#include "TLatex.h"
#include "TColor.h"
#include "TLegend.h"
#include "TLine.h"

#include <iostream>
#include <string>
#include <vector>
#include <map>

// Convert directory name -> human label for header.
// Return "" if not recognized.
static TString CentralityLabelFromDir(const std::string& cname)
{
  if (cname == "CENT_0_10")  return "0-10 %";
  if (cname == "MID_20_40")  return "20-40 %";
  if (cname == "PERI_60_80") return "60-80 %";
  return "";
}

// -------------------- cuts (VISUAL ONLY) --------------------
const double RECO_DUMMY_CUT = -500.0;

// Jet-area cut lines (draw only; do not apply)
const double CUT_AREA_02 = 0.07;  // R = 0.2
const double CUT_AREA_03 = 0.20;  // R = 0.3
const double CUT_AREA_04 = 0.40;  // R = 0.4

// -------------------- utils --------------------
static TString SafeName(const char* name) {
  TString s(name);
  s.ReplaceAll("/", "_");
  s.ReplaceAll(" ", "_");
  return s;
}

static void NormalizeProb(TH1* h) {
  if (!h) return;
  const double I = h->Integral();
  if (I > 0) h->Scale(1.0 / I);
}

static double AreaCutForR(const std::string& rname) {
  if (rname == "R0.2") return CUT_AREA_02;
  if (rname == "R0.3") return CUT_AREA_03;
  if (rname == "R0.4") return CUT_AREA_04;
  return -1.0;
}

// Draw a vertical red dashed line at x = xcut spanning current pad y-range.
// ROOT5/CINT-safe pattern: allocate on heap and do not delete.
static void DrawVLineX(double xcut)
{
  if (!gPad) return;
  if (xcut <= 0) return;

  gPad->Update(); // ensures ranges are computed

  const double y1 = gPad->GetUymin();
  const double y2 = gPad->GetUymax();

  TLine* ln = new TLine(xcut, y1, xcut, y2);
  ln->SetLineColor(kRed);
  ln->SetLineWidth(3);
  ln->SetLineStyle(2);
  ln->Draw("SAME");
}

// -------------------- palette (push to yellow/green, avoid red) --------------------
static void SetYellowGreenPalette() {
  const Int_t NRGBs = 4;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = {0.00, 0.35, 0.70, 1.00};

  Double_t red[NRGBs]   = {0.05, 0.00, 0.20, 1.00};
  Double_t green[NRGBs] = {0.10, 0.80, 0.95, 1.00};
  Double_t blue[NRGBs]  = {0.50, 0.75, 0.25, 0.00};

  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
}

// -------------------- style --------------------
static void BasicStyle() {
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetTitle(0);

  // remove ticks on top/right
  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);
  gStyle->SetErrorX(0);        // remove horizontal error bars on TH1
  gStyle->SetEndErrorSize(2);  // smaller caps (E1). Still useful elsewhere.
  SetYellowGreenPalette();
}

// -------------------- header text --------------------
enum EHeaderCorner { kTopRight = 0, kTopLeft = 1 };

static void DrawHeader(EHeaderCorner corner = kTopRight,
                       const TString& centLabel = "")
{
  if (!gPad) return;

  const double yTop = 1.0 - gPad->GetTopMargin() - 0.02;
  const double dy   = 0.050;

  double x = 0.0;
  int align = 13; // left-top

  if (corner == kTopRight) {
    x = 1.0 - gPad->GetRightMargin() - 0.02;
    align = 33; // right-top
  } else {
    x = gPad->GetLeftMargin() + 0.02;
    align = 13; // left-top
  }

  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(align);

  tx.SetTextSize(0.040);

  if (centLabel.Length() > 0) {
    tx.DrawLatex(x, yTop,
      Form("STAR Au+Au %s  #sqrt{#it{s}_{NN}} = 200 GeV", centLabel.Data()));
  } else {
    tx.DrawLatex(x, yTop,
      "STAR Au+Au  #sqrt{#it{s}_{NN}} = 200 GeV");
  }

  tx.SetTextSize(0.034);
  tx.DrawLatex(x, yTop - dy, "THIS THESIS");
}

// -------------------- axis cosmetics helpers --------------------
static void ApplyAxisStyle1D(TH1* h)
{
  if (!h) return;

  // modest increase vs ROOT defaults
  const double lab = 0.045;  // label size
  const double tit = 0.050;  // title size

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

static void ApplyAxisStyle2D(TH2* h)
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
  if (h->GetZaxis()) {
    h->GetZaxis()->SetLabelSize(lab);
    h->GetZaxis()->SetTitleSize(tit);
    h->GetZaxis()->SetTitleOffset(1.10);
  }
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


// -------------------- saving helpers --------------------
static void Save1D(TH1* hist,
                   const TString& outPdf,
                   const char* xTitle,
                   const char* yTitle,
                   bool logY,
                   const char* drawOpt,
                   EHeaderCorner headerCorner,
                   double vline_x = -1.0,
                   const TString& centLabel = "",
                   int cw = 800, int ch = 600)

{
  if (!hist) return;

  TH1* h = (TH1*)hist->Clone(Form("%s_plot", hist->GetName()));
  if (!h) return;
  h->SetDirectory(0);
  h->SetTitle("");

  if (h->GetXaxis()) h->GetXaxis()->SetTitle(xTitle ? xTitle : "");
  if (h->GetYaxis()) h->GetYaxis()->SetTitle(yTitle ? yTitle : "");

  ApplyAxisStyle1D(h);

  TCanvas c(Form("c_%s", h->GetName()), "", cw, ch);

  // increased margins to accommodate larger labels/titles
  c.SetLeftMargin(0.15);
  c.SetBottomMargin(0.15);
  c.SetTopMargin(0.06);
  c.SetRightMargin(0.05);

  c.SetTickx(0);
  c.SetTicky(0);

  if (logY) {
    c.SetLogy();
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.9);
    h->SetLineWidth(1);
  }

  if (logY) {
    const double minPos = SmallestPositiveBin(h);
    if (minPos > 0) h->SetMinimum(0.5 * minPos);  // or 0.8*minPos
  }

  TString opt = (drawOpt && drawOpt[0]) ? drawOpt : "E1";
  if (logY && opt.Contains("E")) {
    // E0 = no end-caps; add P so you see points, not just bars
    opt = "E0 P";
  }
  h->Draw(opt.Data());

  c.Update();                  // important on log pads
  if (vline_x > 0) DrawVLineX(vline_x);
  DrawHeader(headerCorner, centLabel);

  c.SaveAs(outPdf.Data());
  delete h;
}

static void Save2D(TH2* hist,
                   const TString& outPdf,
                   const char* xTitle,
                   const char* yTitle,
                   bool logZ,
                   EHeaderCorner headerCorner,
                   double vline_x = -1.0,
                   const TString& centLabel = "",
                   int cw = 800, int ch = 600)

{
  if (!hist) return;

  SetYellowGreenPalette();

  TH2* h = (TH2*)hist->Clone(Form("%s_plot", hist->GetName()));
  if (!h) return;
  h->SetDirectory(0);
  h->SetTitle("");

  if (h->GetXaxis()) h->GetXaxis()->SetTitle(xTitle ? xTitle : "");
  if (h->GetYaxis()) h->GetYaxis()->SetTitle(yTitle ? yTitle : "");

  ApplyAxisStyle2D(h);

  TCanvas c(Form("c_%s", h->GetName()), "", cw, ch);

  // increased left/bottom; keep right larger for colorbar
  c.SetLeftMargin(0.15);
  c.SetBottomMargin(0.15);
  c.SetTopMargin(0.06);
  c.SetRightMargin(0.16); // colorbar

  c.SetTickx(0);
  c.SetTicky(0);

  if (logZ) c.SetLogz();

  h->Draw("COLZ");
  c.Update();
  if (vline_x > 0) DrawVLineX(vline_x);
  DrawHeader(headerCorner, centLabel);

  c.SaveAs(outPdf.Data());
  delete h;
}

// -------------------- main --------------------
void make_qa_pdfs(const char *infile = "data_merged.root")
{
  BasicStyle();

  gSystem->mkdir("pdf/QA", kTRUE);
  gSystem->mkdir("pdf/Jets", kTRUE);

  TFile* fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Error: cannot open input file " << infile << "\n";
    return;
  }

  std::cout << "Input : " << infile << "\n";

  // ============================================================
  // 1) QA histograms (TList named "QA_histograms")
  // ============================================================
  {
    TObject* qaObj = fin->Get("QA_histograms");
    TList* qaList = qaObj ? dynamic_cast<TList*>(qaObj) : 0;

    if (!qaList) {
      std::cerr << "QA_histograms not found as TList. Got: "
                << (qaObj ? qaObj->ClassName() : "null") << "\n";
    } else {
      TH1* hcent9   = (TH1*)qaList->FindObject("hcent9");
      TH1* hzvertex = (TH1*)qaList->FindObject("hzvertex");
      TH1* hTrackPt = (TH1*)qaList->FindObject("hTrackPt");
      TH1* hTowE    = (TH1*)qaList->FindObject("hTowE");
      TH1* hDca     = (TH1*)qaList->FindObject("hTrackDca");
      TH1* hEvt0    = (TH1*)qaList->FindObject("hEventStat0");
      TH1* hEvt1    = (TH1*)qaList->FindObject("hEventStat1");

      TH2* hTowEtaPhi   = (TH2*)qaList->FindObject("hTowEtaPhi");
      TH2* hTrackEtaPhi = (TH2*)qaList->FindObject("hTrackEtaPhi");
      TH2* hRhoVsRefMult = (TH2*)qaList->FindObject("hRhoVsRefMult");

      // --- Requested axis labels ---
      if (hcent9) Save1D(hcent9, "pdf/QA/hcent9.pdf",
                         "Centrality", "Counts", false, "HIST", kTopRight);

      // EventStat: add y-axis label = Counts (x left empty by request)
      if (hEvt0)  Save1D(hEvt0,  "pdf/QA/hEventStat0.pdf",
                         "", "Counts", false, "HIST", kTopRight);
      if (hEvt1)  Save1D(hEvt1,  "pdf/QA/hEventStat1.pdf",
                         "", "Counts", false, "HIST", kTopRight);

      if (hzvertex) Save1D(hzvertex, "pdf/QA/hzvertex.pdf",
                           "#it{v}_{z} (cm)", "Weighted counts", false, "E1", kTopRight);
      if (hTrackPt) Save1D(hTrackPt, "pdf/QA/hTrackPt.pdf",
                           "#it{p}_{T}^{track} (GeV/#it{c})", "Weighted counts", true, "E1", kTopRight);
      if (hTowE)    Save1D(hTowE, "pdf/QA/hTowE.pdf",
                           "#it{E}_{T} (GeV)", "Weighted counts", true, "E1", kTopRight);
      if (hDca)     Save1D(hDca, "pdf/QA/hTrackDca.pdf",
                           "DCA (cm)", "Counts", true, "E1", kTopRight);

      if (hTowEtaPhi)   Save2D(hTowEtaPhi,   "pdf/QA/hTowEtaPhi.pdf",   "#eta_{tow} (-)", "#varphi_{tow} (rad)", true, kTopRight);
      if (hTrackEtaPhi) Save2D(hTrackEtaPhi, "pdf/QA/hTrackEtaPhi.pdf", "#eta_{trk} (-)", "#varphi_{trk} (rad)", true, kTopRight);
      if (hRhoVsRefMult) {
        Save2D(hRhoVsRefMult, "pdf/QA/hRhoVsRefMult.pdf",
               "refMult", "#rho (GeV/#it{c} per unit area)", true, kTopRight);
      }
    }
  }

  // ============================================================
  // 2) Jets: loop R*/CENT*/JetTree
  // ============================================================
  TIter itR(fin->GetListOfKeys());
  TKey* keyR = 0;

  while ((keyR = (TKey*)itR())) {
    if (std::string(keyR->GetClassName()) != "TDirectoryFile") continue;

    std::string rname = keyR->GetName();
    if (rname.empty() || rname[0] != 'R') continue;

    TDirectoryFile* dirR = (TDirectoryFile*)fin->Get(rname.c_str());
    if (!dirR) continue;

    TString pdfR = Form("pdf/Jets/%s", rname.c_str());
    gSystem->mkdir(pdfR.Data(), kTRUE);

    const double areaCutLine = AreaCutForR(rname);

    TIter itC(dirR->GetListOfKeys());
    TKey* keyC = 0;

    while ((keyC = (TKey*)itC())) {
      if (std::string(keyC->GetClassName()) != "TDirectoryFile") continue;

      std::string cname = keyC->GetName();
      TDirectoryFile* dirC = (TDirectoryFile*)dirR->Get(cname.c_str());
      if (!dirC) continue;

      TString centLabel = CentralityLabelFromDir(cname);

      TTree* tree = (TTree*)dirC->Get("JetTree");
      if (!tree) continue;

      TString pdfC = Form("%s/%s", pdfR.Data(), cname.c_str());
      gSystem->mkdir(pdfC.Data(), kTRUE);

      // selection (NO area cut applied here)
      TString sel = Form("(reco_pt > %g) && (reco_trigger_match==1)", RECO_DUMMY_CUT);
      TString cut = Form("(centralityWeight) * (%s)", sel.Data());

      TString tag = Form("%s_%s", rname.c_str(), cname.c_str());

      TH1D* hAreaJet = new TH1D(Form("hArea_%s", tag.Data()), "", 120, 0.0, 1.0);
      TH1D* hNEFJet  = new TH1D(Form("hNEF_%s",  tag.Data()), "", 120, 0.0, 1.0);
      TH2D* hPtArea  = new TH2D(Form("h2_ptcorr_vs_area_%s", tag.Data()),
                                "", 120, 0.0, 1.0, 140, -20.0, 80.0);

      hAreaJet->Sumw2();
      hNEFJet->Sumw2();
      hPtArea->Sumw2();

      tree->Draw(Form("reco_area>>%s", hAreaJet->GetName()), cut.Data(), "goff");
      tree->Draw(Form("reco_neutral_fraction>>%s", hNEFJet->GetName()), cut.Data(), "goff");
      tree->Draw(Form("reco_pt_corr:reco_area>>%s", hPtArea->GetName()), cut.Data(), "goff");

      NormalizeProb(hNEFJet);

      // 1D area: draw red cut line at A = CUT_AREA_R
      Save1D(hAreaJet,
             Form("%s/%s.pdf", pdfC.Data(), SafeName(hAreaJet->GetName()).Data()),
             "A (-)", "Weighted entries", true, "E1", kTopRight,
             areaCutLine, centLabel);

      // 1D NEF: unchanged
      Save1D(hNEFJet,
             Form("%s/%s.pdf", pdfC.Data(), SafeName(hNEFJet->GetName()).Data()),
             "NEF (-)", "Probability (a.u.)", true, "E1", kTopLeft, -1.0, centLabel);

      // 2D ptcorr vs area: draw red cut line at same A
      Save2D(hPtArea,
             Form("%s/%s.pdf", pdfC.Data(), SafeName(hPtArea->GetName()).Data()),
             "A (-)", "#it{p}_{T,reco}^{corr} (GeV/#it{c})", true, kTopRight,
             areaCutLine, centLabel);

      delete hAreaJet;
      delete hNEFJet;
      delete hPtArea;
    }
  }

  fin->Close();
  std::cout << "Done. PDFs under pdf/QA, pdf/Jets\n";
}
