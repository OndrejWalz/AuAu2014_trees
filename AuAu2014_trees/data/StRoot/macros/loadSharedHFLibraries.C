#include "TSystem.h"

void loadSharedHFLibraries() {

  gROOT->LoadMacro(
      "$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  // TString fastjet =
  //     "/cvmfs/star.sdcc.bnl.gov/star-spack/spack/opt/spack/linux-rhel7-x86_64/"
  //     "gcc-4.8.5/fastjet-3.3.4-j5evuymea6juu4tkqxxim6nj3z6ldbg3";
  TString ld = gSystem->Getenv("LD_LIBRARY_PATH");
  if (!ld.Contains("fastjet")) {
    cout << "LD_LIBRARY_PATH does not contain fastjet, please set it up"
         << endl;
    return;
  }

  gSystem->Load("libfastjet");
  gSystem->Load("libfastjetplugins");
  gSystem->Load("libfastjettools");

  gSystem->Load("StBTofUtil");
  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StPicoCuts");
  gSystem->Load("StPicoJetMaker");
  gSystem->Load("StRefMultCorr");
  gSystem->Load("St_db_Maker");
  gSystem->Load("StDaqLib");
  gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcADCtoEMaker");
  gSystem->Load("StEpcMaker");
  gSystem->Load("StTriggerUtilities");
  gSystem->Load("StDbBroker");
  gSystem->Load("libgeometry_Tables"); // rember, order of loading makers
                                       // matters

  cout << " loading of shared HF libraries are done" << endl;

  // -->>> ADD your own library/class HERE
  gSystem->Load("StPicoHFJetMaker"); // analysis maker should be last
}
