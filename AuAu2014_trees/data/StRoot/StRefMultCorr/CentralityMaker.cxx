//----------------------------------------------------------------------------------------------------
// $Id: CentralityMaker.cxx,v 1.5 2015/05/22 06:51:56 hmasui Exp $
// $Log: CentralityMaker.cxx,v $
//
//
//----------------------------------------------------------------------------------------------------

#include <iostream>
#include "StRefMultCorr.h"
#include "CentralityMaker.h"

using std::cout ;
using std::endl ;

ClassImp(CentralityMaker)

  CentralityMaker* CentralityMaker::fInstance = 0 ;

//____________________________________________________________________________________________________
CentralityMaker::CentralityMaker()
{
  // Create instance for centrality classes
  fRefMultCorr_P18ih_VpdMB30_MidLow = new StRefMultCorr("refmult_P18ih_VpdMB30_MidLow") ;
}

//____________________________________________________________________________________________________
CentralityMaker::~CentralityMaker()
{ }

//____________________________________________________________________________________________________
CentralityMaker* CentralityMaker::instance()
{
  if ( !fInstance ) {
    // Initialize StRefMultCorr only once
    fInstance = new CentralityMaker() ;
  }

  return fInstance ;
}

//____________________________________________________________________________________________________
StRefMultCorr* CentralityMaker::getRefMultCorr_P18ih_VpdMB30_MidLow()
{
    return fRefMultCorr_P18ih_VpdMB30_MidLow;
}

//____________________________________________________________________________________________________
void CentralityMaker::help() const
{
  cout << endl;
  cout << "//------------------------------------------------------------------------------" << endl;
  cout << "How to get centrality bins by CentralityMaker ?" << endl;
  cout << "  (Please also take a look at StRoot/StRefMultCorr/macros/getCentralityBins.C" << endl;
  cout << endl;
  cout << "1. Initialize run index to read proper data base" << endl;
  cout << "  Suppose we want to get centrality from refmult at run index 11078000" << endl;
  cout << endl;
  cout << "  // NOTE:" << endl;
  cout << "  //  Use BBC coincidence rate (NOT ZDC coincidence rate) for refmult2)" << endl;
  cout << "  StRefMultCorr* refmultCorr = CentralityMaker::instance()->getRefMultCorr();" << endl;
  cout << "  refmultCorr->init(11078000);" << endl;
  cout << endl;
  cout << "2. Initialize relevant variables event-by-event" << endl;
  cout << endl;
  cout << "  // NOTE:" << endl;
  cout << "  //  1st argument is original multiplicity" << endl;
  cout << "  //  If one wants to have centrality from refmult2, you have to put refmult2" << endl;
  cout << "  refmultCorr->initEvent(refmult, vz, zdcCoincidenceRate);" << endl;
  cout << endl;
  cout << "3. Get centrality bins" << endl;
  cout << endl;
  cout << "  const Int_t cent16 = refmultCorr->getCentralityBin16() ;" << endl;
  cout << "  const Int_t cent9  = refmultCorr->getCentralityBin9() ;" << endl;
  cout << "//------------------------------------------------------------------------------" << endl;
  cout << endl;
}


