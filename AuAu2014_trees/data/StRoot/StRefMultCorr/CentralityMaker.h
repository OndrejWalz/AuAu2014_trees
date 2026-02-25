//----------------------------------------------------------------------------------------------------
//  * Interface of StRefMultCorr for possible extention of StRefMultCorr class to the other 
//    centrality measure, such as refmult2.
//  * This interface is also useful when StRefMultCorr needs to be called from two or more different 
//    makers in order to have exactly the same corrected refmult and centrality bins among different makers.
//
//  There is only one change you have to make
//    Replace
//      StRefMultCorr* refmultCorr = new StRefMultCorr();
//    to
//      StRefMultCorr* refmultCorr = CentralityMaker::instance()->getRefMultCorr();
//
//  authors: Hiroshi Masui
//----------------------------------------------------------------------------------------------------

#ifndef __CentralityMaker_h__
#define __CentralityMaker_h__

class StRefMultCorr ;
#include "Rtypes.h"

//____________________________________________________________________________________________________
class CentralityMaker {
  public:
    static CentralityMaker* instance(); // Use this function to access StRefMultCorr
    virtual ~CentralityMaker(); /// Default destructor

    // Interface
    StRefMultCorr* getRefMultCorr_P18ih_VpdMB30_MidLow()  ; // for P18ih, VPDMB-30; |vz| < 30, Mid and low luminosity

    // Print help messages
    void help() const ;

  private:
    CentralityMaker() ; // Constructor is private
    static CentralityMaker* fInstance ; // Static pointer of CentralityMaker

    // Centrality correction classes
    StRefMultCorr* fRefMultCorr_P18ih_VpdMB30_MidLow; // for P18ih, VPDMB-30; |vz| < 30,  Mid and low luminosity
 
    ClassDef(CentralityMaker, 0)
};
#endif

