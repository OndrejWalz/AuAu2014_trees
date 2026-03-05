#include "TFile.h"
#include "TDirectoryFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TH1I.h"
#include "TString.h"

#include <map>
#include <vector>
#include <iostream>
#include <math.h>  // for fabs

struct EventInfo {
  int runId;
  int eventId;
};

// helper: quantize mc_sum_pt to integer key (avoid float epsilon)
// 0.001 GeV resolution -> multiply by 1000 and round manually
long long make_key(float sum_pt) {
  double x = (double)sum_pt * 1000.0; // 0.001 GeV resolution
  if (x >= 0.0)
    return (long long)(x + 0.5);
  else
    return (long long)(x - 0.5);
}

void check_mc_sum_pt(const char* infile = "embedding_merged.root")
{
  TFile* f = TFile::Open(infile, "READ");
  if (!f || f->IsZombie()) {
    std::cerr << "Cannot open " << infile << std::endl;
    return;
  }

  std::cout << "Checking mc_sum_pt duplicates in " << infile << std::endl;

  // --- Find the first JetTree (any R / centrality) ---
  TTree* tree = 0;

  TIter nextR(f->GetListOfKeys());
  TKey* keyR = 0;
  while ((keyR = (TKey*)nextR())) {
    if (std::string(keyR->GetClassName()) != "TDirectoryFile") continue;
    TDirectoryFile* dirR = (TDirectoryFile*)keyR->ReadObj();
    if (!dirR) continue;

    TIter nextC(dirR->GetListOfKeys());
    TKey* keyC = 0;
    while ((keyC = (TKey*)nextC())) {
      if (std::string(keyC->GetClassName()) != "TDirectoryFile") continue;
      TDirectoryFile* dirC = (TDirectoryFile*)keyC->ReadObj();
      if (!dirC) continue;

      tree = (TTree*)dirC->Get("JetTree");
      if (tree) {
        std::cout << "Using tree " << dirR->GetName()
                  << "/" << dirC->GetName() << "/JetTree" << std::endl;
        break;
      }
    }
    if (tree) break;
  }

  if (!tree) {
    std::cerr << "No JetTree found!" << std::endl;
    f->Close();
    return;
  }

  // --- Set up branches ---
  int   runId     = 0;
  int   eventId   = 0;
  float mc_sum_pt = 0.0f;

  if (!tree->GetBranch("runId") ||
      !tree->GetBranch("eventId") ||
      !tree->GetBranch("mc_sum_pt")) {
    std::cerr << "Tree is missing runId/eventId/mc_sum_pt branches." << std::endl;
    f->Close();
    return;
  }

  tree->SetBranchAddress("runId",     &runId);
  tree->SetBranchAddress("eventId",   &eventId);
  tree->SetBranchAddress("mc_sum_pt", &mc_sum_pt);

  // key -> count, and key -> list of events
  std::map<long long,int> counts;
  std::map<long long, std::vector<EventInfo> > eventsForKey;

  // ensure we only count each (runId,eventId) once = per event
  std::map< std::pair<int,int>, bool > seenEvent;

  Long64_t nentries = tree->GetEntries();
  std::cout << "Tree entries: " << nentries << std::endl;

  Long64_t i;
  for (i = 0; i < nentries; ++i) {
    tree->GetEntry(i);

    std::pair<int,int> evKey(runId, eventId);
    if (seenEvent.find(evKey) != seenEvent.end()) {
      // already processed this event
      continue;
    }
    seenEvent[evKey] = true;

    long long key = make_key(mc_sum_pt);
    counts[key]++;

    EventInfo info;
    info.runId   = runId;
    info.eventId = eventId;
    eventsForKey[key].push_back(info);
  }

  std::cout << "Unique events counted: " << seenEvent.size() << std::endl;
  std::cout << "Unique mc_sum_pt keys: " << counts.size() << std::endl;

  // --- Build multiplicity histogram: how many sums appear N times? ---
  int maxMult = 0;
  std::map<long long,int>::const_iterator it;
  for (it = counts.begin(); it != counts.end(); ++it) {
    if (it->second > maxMult) maxMult = it->second;
  }

  TH1I* hMultiplicity =
    new TH1I("hMultiplicity",
             "Multiplicity of identical mc_sum_pt;"
             "number of events with same mc_sum_pt;number of distinct sums",
             maxMult+1, -0.5, maxMult+0.5);

  int nDuplicates = 0;
  std::map<long long,int>::const_iterator it2;
  for (it2 = counts.begin(); it2 != counts.end(); ++it2) {
    int mult = it2->second;
    hMultiplicity->Fill(mult);
    if (mult > 1) nDuplicates++;
  }

  std::cout << "Number of distinct mc_sum_pt values that repeat (mult>1): "
            << nDuplicates << std::endl;

  // Print a few example duplicates
  int printed = 0;
  std::map<long long,int>::const_iterator it3;
  for (it3 = counts.begin(); it3 != counts.end(); ++it3) {
    if (it3->second <= 1) continue;

    double sumPt = ((double)it3->first) / 1000.0;
    std::cout << "mc_sum_pt = " << sumPt << " GeV appears "
              << it3->second << " times" << std::endl;

    std::vector<EventInfo>& evs = eventsForKey[it3->first];
    for (size_t ie = 0; ie < evs.size(); ++ie) {
      std::cout << "   run " << evs[ie].runId
                << ", event " << evs[ie].eventId << std::endl;
    }

    printed++;
    if (printed >= 10) break;
  }

  TFile fout("mc_sum_pt_check.root", "RECREATE");
  hMultiplicity->Write();
  fout.Close();

  std::cout << "Wrote multiplicity histogram to mc_sum_pt_check.root" << std::endl;

  f->Close();
}