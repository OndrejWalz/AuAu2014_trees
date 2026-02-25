// make_response_matrix_pdfs.C (ROOT5/6, STAR-friendly)

#include "TFile.h"
#include "TKey.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TH1.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TString.h"
#include "TROOT.h"
#include "TCollection.h"

#include <iostream>

// ---------------- style ----------------
static void ApplyNiceStyle()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  gStyle->SetTitleSize(0.045, "XYZ");
  gStyle->SetLabelSize(0.040, "XYZ");
  gStyle->SetTitleOffset(1.15, "X");
  gStyle->SetTitleOffset(1.20, "Y");
  gStyle->SetTitleOffset(1.10, "Z");

  gStyle->SetNumberContours(255);
}

// ---------------- helpers ----------------
static TString NiceCentLabel(const TString& centTag)
{
  if (centTag == "CENT_0_10")   return "0-10 %";
  if (centTag == "CENT_10_20")  return "10-20 %";
  if (centTag == "CENT_20_40")  return "20-40 %";
  if (centTag == "CENT_40_60")  return "40-60 %";
  if (centTag == "CENT_60_80")  return "60-80 %";

  if (centTag == "MID_20_40")   return "20-40 %";
  if (centTag == "PERI_60_80")  return "60-80 %";

  // fallback: parse *_X_Y
  Ssiz_t u1 = centTag.Index("_");
  if (u1 != kNPOS) {
    TString rest = centTag(u1+1, centTag.Length()-u1-1);
    Ssiz_t u2 = rest.Index("_");
    if (u2 != kNPOS) {
      TString a = rest(0,u2);
      TString b = rest(u2+1, rest.Length()-u2-1);
      if (a.IsDigit() && b.IsDigit()) return Form("%s--%s %%", a.Data(), b.Data());
    }
  }
  return centTag;
}

static void ParseFolderTags(const TString& dirName,
                            double& Rval, double& ptLeadVal,
                            TString& centTag)
{
  Rval = -1.0;
  ptLeadVal = -1.0;
  centTag = "";

  // R0.x
  Ssiz_t pR = dirName.Index("R0.");
  if (pR != kNPOS) {
    TString tmp = dirName(pR, 4); // "R0.2"
    TString rnum = tmp; rnum.ReplaceAll("R","");
    Rval = rnum.Atof();
  }

  // centrality
  if (dirName.Contains("CENT_0_10")) centTag="CENT_0_10";
  else if (dirName.Contains("CENT_10_20")) centTag="CENT_10_20";
  else if (dirName.Contains("CENT_20_40")) centTag="CENT_20_40";
  else if (dirName.Contains("CENT_40_60")) centTag="CENT_40_60";
  else if (dirName.Contains("CENT_60_80")) centTag="CENT_60_80";
  else if (dirName.Contains("MID_20_40")) centTag="MID_20_40";
  else if (dirName.Contains("PERI_60_80")) centTag="PERI_60_80";

  // ptleadX
  Ssiz_t ppl = dirName.Index("ptlead");
  if (ppl != kNPOS) {
    TString rest = dirName(ppl, dirName.Length()-ppl);
    rest.ReplaceAll(";1","");
    Ssiz_t u = rest.Index("_");
    if (u != kNPOS) rest = rest(0,u); // "ptlead5"
    TString plnum = rest; plnum.ReplaceAll("ptlead","");
    ptLeadVal = plnum.Atof();
  }
}

static TH2* FindResponseTH2(TDirectory* d, const TString& preferredName)
{
  if (!d) return 0;

  TObject* obj = d->Get(preferredName);
  if (obj && obj->InheritsFrom(TH2::Class())) return (TH2*)obj;

  const char* fallbacks[] = {
    "hRespRecoVsTruth_full",
    "hResponseRecoVsTruth_full",
    "hRespRecoVsTruth"
  };
  for (unsigned i=0;i<sizeof(fallbacks)/sizeof(fallbacks[0]);++i) {
    obj = d->Get(fallbacks[i]);
    if (obj && obj->InheritsFrom(TH2::Class())) return (TH2*)obj;
  }

  // first TH2 key
  TIter nextkey(d->GetListOfKeys());
  TKey* key = 0;
  while ((key = (TKey*)nextkey())) {
    TClass* cl = gROOT->GetClass(key->GetClassName());
    if (!cl) continue;
    if (cl->InheritsFrom(TH2::Class())) {
      TObject* h = d->Get(key->GetName());
      if (h && h->InheritsFrom(TH2::Class())) return (TH2*)h;
    }
  }
  return 0;
}

// ---------------- plotting ----------------
static void DrawOne(TH2* h2,
                    const TString& outPdf,
                    const TString& dname,
                    double Rval, double ptLeadVal,
                    const TString& centTag,
                    double xmin, double xmax,
                    double ymin, double ymax,
                    int idx)
{
  if (!h2) return;

  // unique clone name is important in ROOT5 batch loops
  TString cname = Form("h_tmp_resp_%d", idx);
  TH2* h = (TH2*)h2->Clone(cname);
  h->SetDirectory(0);

  h->GetXaxis()->SetTitle("#it{p}_{T}^{meas} (GeV/c)");
  h->GetYaxis()->SetTitle("#it{p}_{T}^{true} (GeV/c)");
  h->GetZaxis()->SetTitle("Probability");

  if (xmax > xmin) h->GetXaxis()->SetRangeUser(xmin, xmax);
  if (ymax > ymin) h->GetYaxis()->SetRangeUser(ymin, ymax);

  h->SetMinimum(1e-6);

  TString canvName = Form("c_resp_%d", idx);
  TCanvas* c = new TCanvas(canvName, canvName, 900, 800);

  c->SetLeftMargin(0.12);
  c->SetBottomMargin(0.12);
  c->SetTopMargin(0.07);
  c->SetRightMargin(0.14);
  c->SetLogz(1);

  h->Draw("COLZ");

  // Title at top-left (not inside the matrix)
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.045);
  lat.DrawLatex(0.12, 0.97, "");

  // --- transparent text box in the left empty x-range region ---
  // This assumes you have xmin < 0 so there is empty space on the left.
  TPaveText* p = new TPaveText(0.14, 0.14, 0.44, 0.44, "NDC"); // left remember
  p->SetFillStyle(0);     // transparent
  p->SetBorderSize(0);
  p->SetTextFont(42);
  p->SetTextSize(0.032);
  p->SetTextAlign(12);

  p->AddText("Response Matrix");
  p->AddText("THIS THESIS");
  p->AddText(Form("anti-#it{k}_{T}, #it{R} = %.1f", (Rval>0?Rval:0.0)));
  p->AddText(Form("#it{p}_{T}^{lead} #geq %.0f GeV/c", (ptLeadVal>=0?ptLeadVal:0.0)));
  p->AddText("#sqrt{s_{NN}} = 200 GeV");
  p->AddText("PYTHIA6 p+p #otimes STAR Au+Au");
  if (centTag.Length()>0)
    p->AddText(Form("%s", NiceCentLabel(centTag).Data()));
  

  p->Draw();

  c->SaveAs(outPdf);

  delete p;
  delete c;
  delete h;
}

// ---------------- recursion ----------------
static void RecurseAndPlot(TDirectory* baseDir,
                           const TString& outDir,
                           const TString& preferredHistName,
                           double xmin, double xmax,
                           double ymin, double ymax,
                           int& idx)
{
  if (!baseDir) return;

  TIter nextkey(baseDir->GetListOfKeys());
  TKey* key = 0;

  while ((key = (TKey*)nextkey())) {
    TClass* cl = gROOT->GetClass(key->GetClassName());
    if (!cl) continue;

    if (cl->InheritsFrom(TDirectory::Class())) {
      TDirectory* d = (TDirectory*) baseDir->Get(key->GetName());
      if (!d) continue;

      TString dname = d->GetName();

      TH2* h2 = FindResponseTH2(d, preferredHistName);
      if (h2) {
        double Rval, ptLeadVal;
        TString centTag;
        ParseFolderTags(dname, Rval, ptLeadVal, centTag);

        TString outPdf = Form("%s/response_%s.pdf", outDir.Data(), dname.Data());
        outPdf.ReplaceAll("/","_"); // safety if any slashes sneak in
        outPdf.ReplaceAll("__","_");

        std::cout << "Plot: dir=" << dname << " -> " << outPdf << std::endl;

        idx++;
        DrawOne(h2, outPdf, dname, Rval, ptLeadVal, centTag, xmin, xmax, ymin, ymax, idx);
      } else {
        RecurseAndPlot(d, outDir, preferredHistName, xmin, xmax, ymin, ymax, idx);
      }
    }
  }
}

void make_response_matrix_pdfs(const char* inFile = "responses_embedding.root",
                               const char* outDir = "pdf/ResponseMatrices",
                               const char* preferredHistName = "hRespRecoVsTruth_full",
                               double xmin = -80.0, double xmax = 60.0,  // wide X by default
                               double ymin = 0.0,    double ymax = 60.0)
{
  ApplyNiceStyle();
  TH1::AddDirectory(kFALSE);
  gROOT->SetBatch(kTRUE);

  gSystem->mkdir(outDir, kTRUE);

  TFile* f = TFile::Open(inFile, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "ERROR: cannot open input file: " << inFile << std::endl;
    return;
  }

  std::cout << "Input:  " << inFile << std::endl;
  std::cout << "Output: " << outDir << std::endl;

  int idx = 0;
  RecurseAndPlot(f, outDir, preferredHistName, xmin, xmax, ymin, ymax, idx);

  f->Close();
  delete f;

  std::cout << "Done." << std::endl;

  // ROOT5/STAR workaround: prevent late teardown segfault in some builds (TPDF unload)
  gSystem->Exit(0);
}
