#include "StPicoCuts.h"

using namespace std;

ClassImp(StPicoCuts)

// _________________________________________________________
StPicoCuts::StPicoCuts() : TNamed("PicoCutsBase", "PicoCutsBase"), 
  mPicoDst(NULL), mEventStatMax(7), mBadRunListFileName("BadRunList_MB.list"), mVzMax(6.), mVzVpdVzMax(3.),
  mNHitsFitMin(20), mRequireHFT(true), mNHitsFitnHitsMax(0.52) {
  
  // -- default constructor

  mPtRange[0] = std::numeric_limits<float>::lowest();
  mPtRange[1] = std::numeric_limits<float>::max();
  mDcaMin = std::numeric_limits<float>::lowest();
  mEta = std::numeric_limits<float>::max();
  mTPCNSigmaMax = std::numeric_limits<float>::max();

}

// _________________________________________________________
StPicoCuts::StPicoCuts(const Char_t *name) : TNamed(name, name), 
  mPicoDst(NULL), mEventStatMax(7), mBadRunListFileName("BadRunList_MB.list"), mVzMax(6.), mVzVpdVzMax(3.),
  mNHitsFitMin(20), mRequireHFT(true), mNHitsFitnHitsMax(0.52) {
  // -- constructor

  mPtRange[0] = std::numeric_limits<float>::lowest();
  mPtRange[1] = std::numeric_limits<float>::max();
  mDcaMin = std::numeric_limits<float>::lowest();
  mEta = std::numeric_limits<float>::max();
  mTPCNSigmaMax = std::numeric_limits<float>::max();

}
// _________________________________________________________
StPicoCuts::~StPicoCuts() { 
  // destructor
}

// _________________________________________________________
void StPicoCuts::initBase() {
  // -- init cuts class

  // -- Read in bad run list and fill vector
  // -----------------------------------------

  // -- open list
  ifstream runs;

  // -- open in working dir
  runs.open(mBadRunListFileName.Data());
  if (!runs.is_open()) {
    runs.open(Form("picoLists/%s", mBadRunListFileName.Data()));
    if (!runs.is_open()) {
      cout << "StPicoCuts::initBase -- Bad run list NOT found :" << mBadRunListFileName << endl;
      cout << "StPicoCuts::initBase -- continue without bad run selection! " << endl;
      //exit(EXIT_FAILURE);
    }
  }

  if (runs.is_open()) {
    Int_t runId = 0;
    while( runs >> runId )
      mVecBadRunList.push_back(runId);
    
    runs.close();

    // -- sort bad runs vector
    std::sort(mVecBadRunList.begin(), mVecBadRunList.end());
  }
	
}

// _________________________________________________________
bool StPicoCuts::isGoodEvent(StPicoDst const * const picoDst, int *aEventCuts) {
  // -- method to check if good event
  //    sets also mPicoDst and mPrimVtx
  
  // -- set current mPicoDst 
  mPicoDst = picoDst;

  // -- get picoDst event
  StPicoEvent* picoEvent = mPicoDst->event();

  // -- set current primary vertex
  mPrimVtx = picoEvent->primaryVertex();

  // -- quick method without providing stats
  if (!aEventCuts) {
		if(mRequireHFT){    
			return (/*isGoodRunHFT(picoEvent) &&*/ isGoodTrigger(picoEvent) &&
		    fabs(picoEvent->primaryVertex().z()) < mVzMax &&
		    fabs(picoEvent->primaryVertex().z() - picoEvent->vzVpd()) < mVzVpdVzMax && picoEvent->refMult() > mRefMultMin && picoEvent->refMult() < mRefMultMax);
		}
		else{
			return (/*isGoodRun(picoEvent) &&*/ isGoodTrigger(picoEvent) &&
		    fabs(picoEvent->primaryVertex().z()) < mVzMax &&
		    fabs(picoEvent->primaryVertex().z() - picoEvent->vzVpd()) < mVzVpdVzMax && picoEvent->refMult() > mRefMultMin && picoEvent->refMult() < mRefMultMax);
		}
  }
    
  // -- reset event cuts
  for (unsigned int ii = 0; ii < mEventStatMax; ++ii)
    aEventCuts[ii] = 0;
  
  unsigned int iCut = 0;

  // -- 0 - before event cuts
  aEventCuts[iCut] = 0;

  // -- 1 - is bad run
  ++iCut;
  /*if(mRequireHFT){
		if(!isGoodRunHFT(picoEvent))
			aEventCuts[iCut] = 1;
	}
	else{*/
		if (!isGoodRun(picoEvent))
    	aEventCuts[iCut] = 1;
	//}

  // -- 2 - No Trigger fired
  ++iCut;
  if (!isGoodTrigger(picoEvent))
    aEventCuts[iCut] = 1;

  // -- 3 - Vertex z outside cut window
  ++iCut;
  if (fabs(picoEvent->primaryVertex().z()) >= mVzMax)
    aEventCuts[iCut] = 1;

  // -- 4 Vertex z - vertex_z(vpd) outside cut window
  ++iCut;
  if (fabs(picoEvent->primaryVertex().z() - picoEvent->vzVpd()) >= mVzVpdVzMax)
    aEventCuts[iCut] = 1;  
  
  // -- 5 check for centrality info
  ++iCut;
  if (picoEvent->refMult() < mRefMultMin || picoEvent->refMult() > mRefMultMax)
      aEventCuts[iCut] = 1;

	/*for (unsigned int ii = 0; ii < mEventStatMax; ++ii)
    cout << aEventCuts[ii] <<"	";
	cout << endl;*/
  // -- is rejected
  bool isGoodEvent = true;
  for (unsigned int ii = 0; ii < mEventStatMax; ++ii)
    if  (aEventCuts[ii])
      isGoodEvent = false;
        
  return isGoodEvent;
}

// _________________________________________________________
bool StPicoCuts::isGoodRun(StPicoEvent const * const picoEvent) const {
  // -- is good run (not in bad runlist)

  return (!(std::binary_search(mVecBadRunList.begin(), mVecBadRunList.end(), picoEvent->runId())));
}

// _________________________________________________________
bool StPicoCuts::isGoodTrigger(StPicoEvent const * const picoEvent) const {
  // -- is good trigger in list of good triggerIds
  for(std::vector<unsigned int>::const_iterator iter = mVecTriggerIdList.begin(); iter != mVecTriggerIdList.end(); ++iter)
    if(picoEvent->isTrigger(*iter))      return true;
	
  return false;
}

// _________________________________________________________
bool StPicoCuts::isGoodTrack(StPicoTrack const * const trk) const {
  // -- require at least one hit on every layer of PXL and (IST or SSD).

    //float nSigma = max(max(fabs(trk->nSigmaPion()), fabs(trk->nSigmaKaon())), fabs(trk->nSigmaProton()));
		float nSigma = min(min(fabs(trk->nSigmaPion()), fabs(trk->nSigmaKaon())), fabs(trk->nSigmaProton()));
    float eta = trk->gMom().PseudoRapidity();
    float dca = (mPrimVtx - trk->origin()).Mag();
		float pT = trk->gMom().Perp();
		float nHitsFit = trk->nHitsFit();

    return ((!mRequireHFT || trk->isHFTTrack()) && nHitsFit > mNHitsFitMin && ((float)nHitsFit/(float)trk->nHitsMax()) > mNHitsFitnHitsMax
            && pT > mPtRange[0] && pT < mPtRange[1] && fabs(eta) < mEta && dca < mDcaMin);

}

// _________________________________________________________
bool StPicoCuts::isGoodPrimaryTrack(StPicoTrack const * const trk) const {

    float eta = trk->pMom().PseudoRapidity();
    float dca = (mPrimVtx - trk->origin()).Mag();
		float pT = trk->pMom().Perp();
		float nHitsFit = trk->nHitsFit();

    return ((!mRequireHFT || trk->isHFTTrack()) && nHitsFit > mNHitsFitMin && ((float)nHitsFit/(float)trk->nHitsMax()) > mNHitsFitnHitsMax
            && pT > mPtRange[0] && pT < mPtRange[1] && fabs(eta) < mEta && dca < mDcaMin);

}

// _________________________________________________________
bool StPicoCuts::isGoodPrimaryTrackWithNsigma(StPicoTrack const * const trk) const {

    //float nSigma = max(max(fabs(trk->nSigmaPion()), fabs(trk->nSigmaKaon())), fabs(trk->nSigmaProton()));
    float nSigma = min(min(fabs(trk->nSigmaPion()), fabs(trk->nSigmaKaon())), fabs(trk->nSigmaProton()));
    float eta = trk->pMom().PseudoRapidity();
    float dca = (mPrimVtx - trk->origin()).Mag();
    float pT = trk->pMom().Perp();
    float nHitsFit = trk->nHitsFit();

    return ((!mRequireHFT || trk->isHFTTrack()) && nHitsFit > mNHitsFitMin && ((float)nHitsFit/(float)trk->nHitsMax()) > mNHitsFitnHitsMax
            && pT > mPtRange[0] && pT < mPtRange[1] && nSigma <= mTPCNSigmaMax && fabs(eta) < mEta && dca < mDcaMin);

}


// _________________________________________________________
bool StPicoCuts::isGoodTowHit(StPicoBTowHit const * const towHit) const {

  	float towE = towHit->energy();
		float towADC = towHit->adc();

    return (towE > mERangeMin && towADC > 4 && towE < mERangeMax);

}

// =======================================================================
