#include "StPicoJetMaker.h"

ClassImp(StPicoJetMaker)

    // _________________________________________________________
    // _________________________________________________________
    StPicoJetMaker::StPicoJetMaker(TString name, StPicoDstMaker *picoMaker,
                                   TString outputBaseFileName)
    : StMaker(name), mPicoDst(NULL), mBField(0.), mOutList(NULL),
      mMakerMode(StPicoJetMaker::kAnalyze), mMcMode(false),
      mOutputTreeName("picoHFtree"), mOutputFileBaseName(outputBaseFileName),
      mPicoDstMaker(picoMaker), mPicoEvent(NULL), mTree(NULL),
      mOutputFileTree(NULL), mOutputFileList(NULL) {
  // -- constructor
}

// _________________________________________________________
StPicoJetMaker::~StPicoJetMaker() {
  // -- destructor

  if (mPicoCuts)
    delete mPicoCuts;
  mPicoCuts = NULL;
}

// _________________________________________________________
Int_t StPicoJetMaker::Init() {
  // -- Inhertited from StMaker
  //    NOT TO BE OVERWRITTEN by daughter class
  //    daughter class should implement InitHF()

  // -- check for cut class
  if (!mPicoCuts)
    mPicoCuts = new StPicoCuts();
  mPicoCuts->init();

  // -- file which holds list of histograms
  mOutputFileList = new TFile(
      Form("%s.%s.root", mOutputFileBaseName.Data(), GetName()), "RECREATE");
  mOutputFileList->SetCompressionLevel(1);

  if (mMakerMode == StPicoJetMaker::kWrite) {
    mOutputFileTree = new TFile(
        Form("%s.%s.root", mOutputFileBaseName.Data(), mOutputTreeName.Data()),
        "RECREATE");
    mOutputFileTree->SetCompressionLevel(1);
    mOutputFileTree->cd();

    // -- create OutputTree
    int BufSize = (int)pow(2., 16.);
    int Split = 1;
    if (!mTree)
      mTree = new TTree("T", "T", BufSize);
    mTree->SetAutoSave(1000000); // autosave every 1 Mbytes
  }

  // -- disable automatic adding of objects to file
  bool oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(false);

  // -- add list which holds all histograms
  mOutList = new TList();
  mOutList->SetName(GetName());
  mOutList->SetOwner(true);

  // -- create event stat histograms
  initializeEventStats();

  // -- call method of daughter class
  InitJets();

  TH1::AddDirectory(oldStatus);

  // -- reset event to be in a defined state
  resetEvent();

  return kStOK;
}

// _________________________________________________________
Int_t StPicoJetMaker::Finish() {
  // -- Inhertited from StMaker
  //    NOT TO BE OVERWRITTEN by daughter class
  //    daughter class should implement FinishHF()

  if (mMakerMode == StPicoJetMaker::kWrite) {
    mOutputFileTree->cd();
    mOutputFileTree->Write();
    mOutputFileTree->Close();
  }
  mOutputFileList->cd();
  mOutList->Write(mOutList->GetName(), TObject::kSingleKey);
  mOutList->Delete();
  // -- call method of daughter class
  FinishJets();

  mOutputFileList->Close();

  return kStOK;
}

// _________________________________________________________
void StPicoJetMaker::resetEvent() {
  // -- reset event
  mIdxPicoParticles.clear();
}

// _________________________________________________________
void StPicoJetMaker::Clear(Option_t *opt) {
  // -- Inhertited from StMaker
  //    NOT TO BE OVERWRITTEN by daughter class
  //    daughter class should implement ClearHF()

  resetEvent();
}

// _________________________________________________________
Int_t StPicoJetMaker::Make() {
  // -- Inhertited from StMaker
  //    NOT TO BE OVERWRITTEN by daughter class
  //    daughter class should implement MakeHF()
  // -- isPion, isKaon, isProton methods are to be
  //    implemented by daughter class (
  //    -> methods of StHFCuts can and should be used

  if (!mPicoDstMaker) {
    LOG_WARN << " StPicoJetMaker - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  mPicoDst = mPicoDstMaker->picoDst();
  if (!mPicoDst) {
    LOG_WARN << " StPicoJetMaker - No PicoDst! Skip! " << endm;
    return kStWarn;
  }
  Int_t iReturn = kStOK;

  if (setupEvent()) {
    UInt_t nTracks = mPicoDst->numberOfTracks();

    // -- Fill vectors of particle types
    if (mMakerMode == StPicoJetMaker::kWrite ||
        mMakerMode == StPicoJetMaker::kAnalyze) {
        TH1D* hDca = (TH1D*)mOutList->FindObject("hTrackDca");
        TH1I* hNHits = (TH1I*)mOutList->FindObject("hTrackNHitsFit");

        for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack) {
          StPicoTrack *trk = mPicoDst->track(iTrack);
          if (!trk) continue;

          double pTtrack = trk->pMom().Perp();
          if (pTtrack > 30) return kStOK;  // keep your veto

          // Fill "full" distributions BEFORE your track cuts
          double dca = (mPrimVtx - trk->origin()).Mag();
          int nHitsFit = trk->nHitsFit();

          if (hDca) hDca->Fill(dca);
          if (hNHits) hNHits->Fill(nHitsFit);

          // Now apply cuts for jet constituents
          if (!mPicoCuts->isGoodPrimaryTrack(trk)) continue;

          mIdxPicoParticles.push_back(iTrack);
        }
    }
    // -- call method of daughter class
    //// Uncomment here to include jets again
    iReturn = MakeJets();

    // TODO: Fill good event histograms here - expand for other observables
    static_cast<TH1I *>(mOutList->FindObject("hevents_acc"))->Fill(1);
    static_cast<TH1I *>(mOutList->FindObject("hrefmult_acc"))
        ->Fill(mPicoDst->event()->refMult());
    static_cast<TH1F *>(mOutList->FindObject("hzvertex_acc"))
        ->Fill(mPrimVtx.z());
    static_cast<TH1D *>(mOutList->FindObject("hdeltaz_acc"))
        ->Fill(mPrimVtx.z() - mPicoDst->event()->vzVpd());
    static_cast<TH2D *>(mOutList->FindObject("hz_refmult_acc"))
        ->Fill(mPrimVtx.z(), mPicoDst->event()->refMult());
    static_cast<TH1I *>(mOutList->FindObject("hrunId_acc"))
        ->Fill(mPicoDst->event()->runId());
  }

  // -- save information about all events, good or bad
  if (mMakerMode == StPicoJetMaker::kWrite)
    mTree->Fill();

  // -- fill basic event histograms - for all events
  // TODO: Fill all event histograms here - expand for other observables
  static_cast<TH1I *>(mOutList->FindObject("hevents"))->Fill(1);
  static_cast<TH1I *>(mOutList->FindObject("hrefmult"))
      ->Fill(mPicoDst->event()->refMult());
  static_cast<TH1F *>(mOutList->FindObject("hzvertex"))->Fill(mPrimVtx.z());
  static_cast<TH1D *>(mOutList->FindObject("hdeltaz"))
      ->Fill(mPrimVtx.z() - mPicoDst->event()->vzVpd());
  static_cast<TH2D *>(mOutList->FindObject("hz_refmult"))
      ->Fill(mPrimVtx.z(), mPicoDst->event()->refMult());
  static_cast<TH1I *>(mOutList->FindObject("hrunId"))
      ->Fill(mPicoDst->event()->runId());

  // -- reset event to be in a defined state
  resetEvent();

  return (kStOK && iReturn);
}

// _________________________________________________________
bool StPicoJetMaker::setupEvent() {

  // -- fill members from pico event, check for good eventa and fill event
  // statistics

  mPicoEvent = mPicoDst->event();

  mBField = mPicoEvent->bField();
  mPrimVtx = mPicoEvent->primaryVertex();

  int aEventStat[mPicoCuts->eventStatMax()];

  bool bResult = mPicoCuts->isGoodEvent(mPicoDst, aEventStat);

  fillEventStats(aEventStat);

  return bResult;
}

// _________________________________________________________
void StPicoJetMaker::initializeEventStats() {
  // -- Initialize event statistics histograms

  int refmultbins = 130;
  float refmultmin = 0;
  float refmultmax = 650;

  int zbins = 250;
  float zmin = -50;
  float zmax = 50;

  const char *aEventCutNames[] = {"all",
                                  "good run",
                                  "trigger",
                                  "#it{v}_{z}",
                                  "#it{v}_{z}-#it{v}^{VPD}_{z}",
                                  "refMult",
                                  "accepted"};

  mOutList->Add(new TH1F(
      "hEventStat0", "Event cut statistics 0;Event Cuts;Events",
      mPicoCuts->eventStatMax(), -0.5, mPicoCuts->eventStatMax() - 0.5));
  TH1F *hEventStat0 = static_cast<TH1F *>(mOutList->Last());

  mOutList->Add(new TH1F(
      "hEventStat1", "Event cut statistics 1;Event Cuts;Events",
      mPicoCuts->eventStatMax(), -0.5, mPicoCuts->eventStatMax() - 0.5));
  TH1F *hEventStat1 = static_cast<TH1F *>(mOutList->Last());

  for (unsigned int ii = 0; ii < mPicoCuts->eventStatMax(); ii++) {
    hEventStat0->GetXaxis()->SetBinLabel(ii + 1, aEventCutNames[ii]);
    hEventStat1->GetXaxis()->SetBinLabel(ii + 1, aEventCutNames[ii]);
  }

  // TODO: Add event ID histograms for AuAu Run 2016 and observables vs day
  // information All event histograms
  mOutList->Add(new TH1I("hevents", "number of events", 2, 0, 2));
  mOutList->Add(new TH1I("hrefmult", "Reference multiplicity", refmultbins,
                         refmultmin, refmultmax));
  mOutList->Add(
      new TH1F("hzvertex", "z-position of primary vertex", zbins, zmin, zmax));
  mOutList->Add(new TH1D("hdeltaz", "zTPC-zVPD; #Delta [cm]", 250, -10, 10));
  mOutList->Add(new TH2D("hz_refmult", "zvertex vs refmult; z [cm]; refMult",
                         zbins, zmin, zmax, refmultbins, refmultmin,
                         refmultmax));

  mOutList->Add(new TH1I("hrunId", "runId", 90913, 15076101, 15167014));

  // Accepted event histograms
  mOutList->Add(new TH1I("hevents_acc", "number of events", 2, 0, 2));
  mOutList->Add(new TH1I("hrefmult_acc", "Reference multiplicity", refmultbins,
                         refmultmin, refmultmax));
  mOutList->Add(new TH1F("hzvertex_acc", "z-position of primary vertex", zbins,
                         zmin, zmax));
  mOutList->Add(new TH1D("hdeltaz_acc", "zTPC-zVPD; #Delta [cm]", 250, -10, 10));
  mOutList->Add(new TH2D("hz_refmult_acc",
                         "zvertex vs refmult; z [cm]; refMult", zbins, zmin,
                         zmax, refmultbins, refmultmin, refmultmax));

  mOutList->Add(new TH1I("hrunId_acc", "accepted events runId", 90913, 15076101,
                         15167014)); // 15076101−15167014
  mOutList->Add(new TH1D("hTrackDca", "Tracks DCA (to PV);DCA (cm);counts", 400, 0.0, 5.0));
  mOutList->Add(new TH1I("hTrackNHitsFit", "Tracks nHitsFit;nHitsFit;counts", 50, 0, 50));

}

//________________________________________________________________________
void StPicoJetMaker::fillEventStats(int *aEventStat) {
  // -- Fill event statistics

  TH1F *hEventStat0 = static_cast<TH1F *>(mOutList->FindObject("hEventStat0"));
  TH1F *hEventStat1 = static_cast<TH1F *>(mOutList->FindObject("hEventStat1"));

  for (unsigned int idx = 0; idx < mPicoCuts->eventStatMax(); ++idx) {
    if (!aEventStat[idx])
      hEventStat0->Fill(idx);
  }

  for (unsigned int idx = 0; idx < mPicoCuts->eventStatMax(); ++idx) {
    if (aEventStat[idx])
      break;
    hEventStat1->Fill(idx);
  }
}
