// usage: root -l -b -q 'print_trig_matched_pairs.C("embedding_pt50_-1.root", 0.2, "CENT_0_10", 50)'


#include "TFile.h"
#include "TDirectoryFile.h"
#include "TTree.h"
#include "TString.h"

#include <iostream>

void print_trig_matched_pairs(const char* infile    = "embedding_merged.root",
                              double      R         = 0.2,
                              const char* centDir   = "CENT_0_10",
                              int         maxPrint  = 200)
{
  // ---------------- open file ----------------
  TFile* f = TFile::Open(infile, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open " << infile << std::endl;
    return;
  }

  // ---------------- locate JetTree for given R, centrality ----------------
  TString rname;
  rname.Form("R%.1f", R);

  TDirectoryFile* dirR = (TDirectoryFile*)f->Get(rname);
  if (!dirR) {
    std::cerr << "Directory " << rname << " not found in file." << std::endl;
    f->Close();
    return;
  }

  TDirectoryFile* dirC = (TDirectoryFile*)dirR->Get(centDir);
  if (!dirC) {
    std::cerr << "Directory " << rname << "/" << centDir << " not found." << std::endl;
    f->Close();
    return;
  }

  TTree* tree = (TTree*)dirC->Get("JetTree");
  if (!tree) {
    std::cerr << "JetTree not found in " << rname << "/" << centDir << std::endl;
    f->Close();
    return;
  }

  std::cout << "Using tree: " << rname << "/" << centDir << "/JetTree" << std::endl;

  // ---------------- set up branches ----------------
  Int_t   runId   = 0;
  Int_t   eventId = 0;

  Float_t mc_pt   = -999.f;
  Float_t mc_eta  = -999.f;
  Float_t mc_phi  = -999.f;

  Float_t reco_pt      = -999.f;
  Float_t reco_pt_corr = -999.f;
  Float_t reco_eta     = -999.f;
  Float_t reco_phi     = -999.f;

  Bool_t  reco_trigger_match = kFALSE;
  Float_t deltaR             = -1.f;

  if (!tree->GetBranch("runId") ||
      !tree->GetBranch("eventId") ||
      !tree->GetBranch("mc_pt")  ||
      !tree->GetBranch("mc_eta") ||
      !tree->GetBranch("mc_phi") ||
      !tree->GetBranch("reco_pt") ||
      !tree->GetBranch("reco_pt_corr") ||
      !tree->GetBranch("reco_eta") ||
      !tree->GetBranch("reco_phi") ||
      !tree->GetBranch("reco_trigger_match") ||
      !tree->GetBranch("deltaR")) {
    std::cerr << "Missing one of the required branches "
              << "(runId,eventId,mc_*,reco_*,reco_trigger_match,deltaR)." << std::endl;
    f->Close();
    return;
  }

  tree->SetBranchAddress("runId",             &runId);
  tree->SetBranchAddress("eventId",           &eventId);
  tree->SetBranchAddress("mc_pt",             &mc_pt);
  tree->SetBranchAddress("mc_eta",            &mc_eta);
  tree->SetBranchAddress("mc_phi",            &mc_phi);
  tree->SetBranchAddress("reco_pt",           &reco_pt);
  tree->SetBranchAddress("reco_pt_corr",      &reco_pt_corr);
  tree->SetBranchAddress("reco_eta",          &reco_eta);
  tree->SetBranchAddress("reco_phi",          &reco_phi);
  tree->SetBranchAddress("reco_trigger_match",&reco_trigger_match);
  tree->SetBranchAddress("deltaR",            &deltaR);

  Long64_t nentries = tree->GetEntries();
  std::cout << "Tree entries: " << nentries << std::endl;

  int nPrinted = 0;
  int nMatchedTrig = 0;

  // ---------------- single pass over tree ----------------
  for (Long64_t i = 0; i < nentries; ++i) {
    tree->GetEntry(i);

    // same logic as in maker: matched jet requires both MC and reco
    Bool_t haveMC   = (mc_pt   >= 0.0f);
    Bool_t haveReco = (reco_pt >= 0.0f);
    Bool_t isMatched     = (haveMC && haveReco);
    Bool_t isTrigMatched = (reco_trigger_match == kTRUE);

    if (!(isMatched && isTrigMatched)) continue;

    nMatchedTrig++;

    if (nPrinted < maxPrint) {
      std::cout << "----------------------------------------\n";
      std::cout << "runId="   << runId
                << "  eventId=" << eventId << "\n";
      std::cout << "  RECO: pt=" << reco_pt
                << "  pt_corr=" << reco_pt_corr
                << "  eta=" << reco_eta
                << "  phi=" << reco_phi
                << "  trigMatch=" << 1 << "\n";
      std::cout << "  MC  : pt=" << mc_pt
                << "  eta=" << mc_eta
                << "  phi=" << mc_phi << "\n";
      std::cout << "  deltaR=" << deltaR << "\n";
      nPrinted++;
    }
  }

  std::cout << "\nTotal jets with (matched && trigger_match) = "
            << nMatchedTrig << std::endl;
  if (nPrinted >= maxPrint) {
    std::cout << "(printout truncated at maxPrint=" << maxPrint << ")\n";
  }

  f->Close();
}