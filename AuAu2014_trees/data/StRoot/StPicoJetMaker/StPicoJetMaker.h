#ifndef StPicoJetMaker_h
#define StPicoJetMaker_h

#include <vector>

#include "StChain/StMaker.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StarClassLibrary/StLorentzVectorF.hh"
#include "TVector3.h"

#include "TChain.h"
#include "TFile.h"
#include "TH2D.h"
#include "TString.h"
#include "TTree.h"

#include "../StPicoCuts/StPicoCuts.h"

class StPicoJetMaker : public StMaker {
public:
  StPicoJetMaker(TString name, StPicoDstMaker *picoMaker,
                 TString outputBaseFileName);
  virtual ~StPicoJetMaker();

  // -- TO BE IMPLEMENTED BY DAUGHTER CLASS
  virtual Int_t InitJets() { return kStOK; }
  virtual Int_t MakeJets() { return kStOK; }
  virtual void ClearJets(Option_t *opt = "") { return; }
  virtual Int_t FinishJets() { return kStOK; }

  void setPicoCuts(StPicoCuts *cuts);
  void setTreeName(const char *tName);
  void setMcMode(bool b);
  void setMakerMode(unsigned short us);

  enum eMakerMode { kAnalyze, kWrite };

  // -- Inhertited from StMaker
  //    NOT TO BE OVERWRITTEN by daughter class
  //    daughter class should implement xxxHF()
  //    -> will be declared as "final" when C++11 is used in STAR
  Int_t Init();
  Int_t Make();
  void Clear(Option_t *opt = "");
  Int_t Finish();

protected:
  bool isMcMode() const;
  unsigned int isMakerMode() const;

  // -- protected members ------------------------

  StPicoDst *mPicoDst;

  StPicoCuts *mPicoCuts;

  float mBField;
  TVector3 mPrimVtx;

  TList *mOutList;

  TTree *mTree; // tree holding output information for further processing

  std::vector<unsigned short> mIdxPicoParticles;

private:
  void resetEvent();
  bool setupEvent();

  void initializeEventStats();
  void fillEventStats(int *aEventStat);

  // -- private members ------------------------

  unsigned int mMakerMode; // use enum of StPicoEventMaker::eMakerMode
  bool mMcMode;            // use MC mode

  TString mOutputTreeName; // name for output trees

  TString
      mOutputFileBaseName; // base name for output files
                           //   for tree     ->
                           //   <mOutputFileBaseName>.<mOutputTreeName>.root for
                           //   histList -> <mOutputFileBaseName>.GetName().root

  StPicoDstMaker *mPicoDstMaker; // ptr to picoDst maker

  StPicoEvent *mPicoEvent; // ptr to picoDstEvent

  TFile *mOutputFileTree; // ptr to file saving the HFtree
  TFile *mOutputFileList; // ptr to file saving the list of histograms
  ClassDef(StPicoJetMaker, 0)
};

inline void StPicoJetMaker::setPicoCuts(StPicoCuts *cuts) { mPicoCuts = cuts; }
inline void StPicoJetMaker::setTreeName(const char *tName) {
  mOutputTreeName = tName;
}

inline void StPicoJetMaker::setMcMode(bool b) { mMcMode = b; }
inline bool StPicoJetMaker::isMcMode() const { return mMcMode; }

inline void StPicoJetMaker::setMakerMode(unsigned short us) { mMakerMode = us; }
inline unsigned int StPicoJetMaker::isMakerMode() const { return mMakerMode; }
#endif
