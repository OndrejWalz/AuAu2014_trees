//  Centrality binning:
//     Bin       Centrality (16)   Centrality (9)
//     0            75-80%            70-80%
//     1            70-75%            60-70%
//     2            65-70%            50-60%
//     3            60-65%            40-50%
//     4            55-60%            30-40%
//     5            50-55%            20-30%
//     6            45-50%            10-20%
//     7            40-45%             5-10%
//     8            35-40%             0- 5%
//     9            30-35%
//    10            25-30%
//    11            20-25%
//    12            15-20%
//    13            10-15%
//    14             5-10%
//    15             0- 5%
//
// ACTUALLY NO, THE 9-BIN CENTRALITY IS ORDERED INVERSELY
//
// Specify the type of multiplicity (default is refmult)
// "refmult"   - reference multiplicity defined in |eta|<0.5
// "refmult2"  - reference multiplicity defined in 0.5<|eta|<1.0
// "refmult3"  - reference multiplicity defined in |eta|<0.5 without protons
// "toftray"   - TOF tray multiplicity
// "grefmult"  - global reference multiplicity defined in |eta|<0.5,dca<3,nHitsFit>10
//------------------------------------------------------------------------------
#ifndef __StRefMultCorr_h__
#define __StRefMultCorr_h__

#include <vector>
#include <map>
#include <iostream>
#include "TString.h"
#include "TRandom.h"

class StRefMultCorr
{
 public:
  StRefMultCorr(const TString name="refmult");
  virtual ~StRefMultCorr();          // Default destructor
  
  void           clear();            // Clear all arrays
  const Char_t*  getTable() const;   // Locate the table text file
  void           read();             // Read input parameters from text file 
  
  void setEvent(int runid, double refmult, double zdc, double vz);    // Give event parameters and calculate centrality
  bool checkEvent(int runid, double refmult, double zdc, double vz);  // Check if the event is okay
  void calculateCentrality(double refmult, double zdc, double vz);    // Calculate centrality

  // Getters
  double refMultCorr() const { return mrefmultcorr; }
  double weight()            { return mweight; }
  int    centrality16()      { return mcentrality_16; }
  int    centrality9()       { return mcentrality_9; }
  double ZDCMin() const      { return mmin_zdc; }
  double ZDCMax() const      { return mmax_zdc; }
  int    runMin() const      { return mmin_run; }
  int    runMax() const      { return mmax_run; }
  double VzMin() const       { return mmin_vz; }
  double VzMax() const       { return mmax_vz; }
  
 private:
  const TString mName ; // refmult, refmult2, refmult3 or toftray (case insensitive)
  
  double mrefmultcorr;
  int mcentrality_16;
  int mcentrality_9;
  double mweight;
  
  int myear;
  int menergy;
  double mmin_vz;
  double mmax_vz;
  double mmin_zdc;
  double mmax_zdc;
  int mmin_run;
  int mmax_run;
  double mweight_bound;
  
  double mvz_norm;
  double mzdc_norm;
  
  std::vector<double> mzdc_par;
  std::vector<double> mvz_par;
  std::vector<double> mweight_par;
  std::vector<unsigned> mcentbin_16;
  std::vector<unsigned> mcentbin_9;
  
  TRandom* gRandom;
  ClassDef(StRefMultCorr, 0)
};
#endif
