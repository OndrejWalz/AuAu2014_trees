// print_event_jets.C
//
// Usage:
//   root -l -b -q 'print_event_jets.C("embedding.root", 0.2, "CENT_0_10", 15076000, 12345)'
//
// Prints all jets in the specified event.

#include "TFile.h"
#include "TDirectoryFile.h"
#include "TTree.h"
#include "TString.h"
#include <iostream>

void print_event_jets(const char* infile   = "embedding_merged.root",
                      double      R        = 0.2,
                      const char* centDir  = "CENT_0_10",
                      int         runReq   = 0,
                      int         eventReq = 0)
{
  // ---------------- open file ----------------
  TFile* f = TFile::Open(infile, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open file: " << infile << std::endl;
    return;
  }

  // ---------------- locate tree ----------------
  TString rname;
  rname.Form("R%.1f", R);

  TDirectoryFile* dirR = (TDirectoryFile*)f->Get(rname);
  if (!dirR) {
    std::cerr << "Directory " << rname << " not found." << std::endl;
    return;
  }

  TDirectoryFile* dirC = (TDirectoryFile*)dirR->Get(centDir);
  if (!dirC) {
    std::cerr << "Directory " << rname << "/" << centDir << " not found." << std::endl;
    return;
  }

  TTree* tree = (TTree*)dirC->Get("JetTree");
  if (!tree) {
    std::cerr << "JetTree not found." << std::endl;
    return;
  }

  std::cout << "Using tree: " << rname << "/" << centDir << "/JetTree\n";

  // ---------------- set up branches ----------------
  Int_t runId = 0, eventId = 0;

  Float_t mc_pt   = -999.f;
  Float_t mc_eta  = -999.f;
  Float_t mc_phi  = -999.f;

  Float_t reco_pt      = -999.f;
  Float_t reco_pt_corr = -999.f;
  Float_t reco_eta     = -999.f;
  Float_t reco_phi     = -999.f;

  Bool_t  reco_trigger_match = kFALSE;
  Float_t deltaR             = -1.f;

  tree->SetBranchAddress("runId", &runId);
  tree->SetBranchAddress("eventId", &eventId);

  if (tree->GetBranch("mc_pt")) tree->SetBranchAddress("mc_pt", &mc_pt);
  if (tree->GetBranch("mc_eta")) tree->SetBranchAddress("mc_eta", &mc_eta);
  if (tree->GetBranch("mc_phi")) tree->SetBranchAddress("mc_phi", &mc_phi);

  tree->SetBranchAddress("reco_pt",      &reco_pt);
  tree->SetBranchAddress("reco_pt_corr", &reco_pt_corr);
  tree->SetBranchAddress("reco_eta",     &reco_eta);
  tree->SetBranchAddress("reco_phi",     &reco_phi);

  tree->SetBranchAddress("reco_trigger_match", &reco_trigger_match);
  tree->SetBranchAddress("deltaR", &deltaR);

  Long64_t nentries = tree->GetEntries();

  std::cout << "Searching for runId=" << runReq
            << " eventId=" << eventReq << " …\n";

  bool found = false;
  int jetIdx = 0;

  // ---------------- loop and print ----------------
  for (Long64_t i = 0; i < nentries; ++i) {
    tree->GetEntry(i);

    if (runId != runReq || eventId != eventReq)
      continue;

    if (!found) {
      std::cout << "\n============================================\n";
      std::cout << "Event found: runId=" << runId
                << " eventId=" << eventId << "\n";
      std::cout << "============================================\n";
      found = true;
    }

    std::cout << "Jet " << jetIdx++ << ":\n";

    // Reco jet info
    std::cout << "  RECO: pt=" << reco_pt
              << "  pt_corr=" << reco_pt_corr
              << "  eta=" << reco_eta
              << "  phi=" << reco_phi
              << "  trigMatch=" << (reco_trigger_match ? 1 : 0)
              << "\n";

    // MC jet info (if available)
    if (mc_pt >= 0) {
      std::cout << "  MC  : pt=" << mc_pt
                << "  eta=" << mc_eta
                << "  phi=" << mc_phi
                << "\n";
    } else {
      std::cout << "  MC  : none\n";
    }

    std::cout << "  deltaR=" << deltaR << "\n";
    std::cout << "--------------------------------------------\n";
  }

  if (!found) {
    std::cout << "No entries found for that run/event.\n";
  }

  f->Close();
}
