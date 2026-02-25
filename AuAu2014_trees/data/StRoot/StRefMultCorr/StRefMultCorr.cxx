#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include "StRefMultCorr.h"
#include "TError.h"
#include "TRandom.h"
#include "TMath.h"

ClassImp(StRefMultCorr)

using namespace std ;

  namespace {
    typedef pair<Double_t, Int_t> keys;
  }

//______________________________________________________________________________
// Default constructor
StRefMultCorr::StRefMultCorr(const TString name)
  : mName(name), mrefmultcorr(-1.0), mcentrality_16(-1), mcentrality_9(-1), mweight(-1.0),
    myear(0), menergy(0), mmin_vz(-999.0), mmax_vz(-999.0), mmin_zdc(0.0), mmax_zdc(60000.0),
    mmin_run(-1), mmax_run(-1), mweight_bound(0.0), mvz_norm(0.0), mzdc_norm(30000.0)
{
  clear();
  
  gRandom = new TRandom();
  read();
}

//______________________________________________________________________________
// Default destructor
StRefMultCorr::~StRefMultCorr()
{
  
}

//______________________________________________________________________________
void StRefMultCorr::clear()
{
  mzdc_par.clear();
  mvz_par.clear();
  mweight_par.clear();
  mcentbin_16.clear();
  mcentbin_9.clear();
}

//______________________________________________________________________________
const Char_t* StRefMultCorr::getTable() const
{
  if ( mName.CompareTo("refmult_P18ih_VpdMB30_MidLow", TString::kIgnoreCase) == 0 ) {
      return "StRoot/StRefMultCorr/Centrality_def_refmult_P18ih_VpdMB30_MidLow.txt";
  }
  else{
    Error("StRefMultCorr::getTable", "No implementation for %s", mName.Data());
    cout << "Current available option is refmult or refmult2 or refmult3 or toftray" << endl;
    return "";
  }
}

//______________________________________________________________________________
void StRefMultCorr::read()
{
  // Open the parameter file and read the data
  const Char_t* inputFileName(getTable());
  ifstream ParamFile(inputFileName);
  if(!ParamFile){
    Error("StRefMultCorr::read", "cannot open %s", inputFileName);
    return;
  }
  cout << "StRefMultCorr::read  Open " << inputFileName << flush ;

  string line ;
  getline(ParamFile,line);

  if(line.find("Start_runId")!=string::npos)
  {
    ParamFile >> myear >> menergy ;
    ParamFile >> mmin_run >> mmax_run >> mmin_vz >> mmax_vz ;
    
    for(Int_t i=0;i<16;i++)
    {
      Int_t centralitybins=-1;
      ParamFile >> centralitybins;
      mcentbin_16.push_back(centralitybins);
      if (i%2==0 || i==15) mcentbin_9.push_back(centralitybins);
    }
    
    ParamFile >> mweight_bound;
  
  for(Int_t i=0;i<7;i++)
  {
    Double_t param=-9999.;
    ParamFile >> param;
    mvz_par.push_back( param );
  }
  for(Int_t i=0;i<7;i++) {
    Double_t param=-9999.;
    ParamFile >> param;
    mweight_par.push_back( param );
  }
  
  for(Int_t i=0;i<2;i++) {
    Double_t param=-9999.;
    ParamFile >> param;
    mzdc_par.push_back( param );
  }
  
  }
 else
  {
    cout << endl;
    Error("StRefMultCorr::read", "Input file is not correct! Wrong structure.");
    return;
  }
  ParamFile.close();

  cout << " [OK]" << endl;
}

//______________________________________________________________________________
void StRefMultCorr::setEvent(int runid, double refmult, double zdc, double vz)
{
  if (checkEvent(runid, refmult, zdc, vz))  calculateCentrality(refmult, zdc, vz);
  else
  {
    cout << "Bad event" << endl;
    mrefmultcorr = refmult;
    mcentrality_9 = -1;
    mcentrality_16 = -1;
    mweight = 0;
  }
}

//______________________________________________________________________________
bool StRefMultCorr::checkEvent(int runid, double refmult, double zdc, double vz)
{
  if (refmult < 0) { cout << "ref mult" << endl; return false; }
  if (runid < mmin_run || runid > mmax_run){ cout << "run" << endl; return false; }
  if (vz < mmin_vz || vz > mmax_vz) { cout << "vz" << endl; return false; }
  if (zdc < mmin_zdc || zdc > mmax_zdc) { cout << "zdc" << endl; return false; }
  return true;
}

//______________________________________________________________________________
void StRefMultCorr::calculateCentrality(double refmult, double zdc, double vz)
{
  // we randomize raw refmult within 1 bin to avoid the peaky structures at low refmult
  double raw_ref = refmult + gRandom->Rndm();
    
  if(mzdc_par.empty() || mvz_par.empty())
  {
    mrefmultcorr = 0.0;
    mcentrality_9 = -1;
    mcentrality_16 = -1;
    mweight = 0.0;
  }
  
  double zdc_scaling = mzdc_par[0] + mzdc_par[1] * zdc / 1000.0;
  double zdc_norm = mzdc_par[0] + mzdc_par[1] * mzdc_norm / 1000.0;
  double zdc_correction = zdc_norm / zdc_scaling;

  double vz_scaling = 0.0;
  double vz_norm = 0.0;
  
  for (int i = 0; i < 7; ++i)
  {
    vz_scaling += mvz_par[i] * pow(vz, i);
    vz_norm += mvz_par[i] * pow(mvz_norm, i);
  }

  double vz_correction = 1.0;
  if (vz_scaling > 0.0)
  {
    vz_correction = vz_norm / vz_scaling;
  }

  mrefmultcorr = raw_ref * vz_correction * zdc_correction;
  
  // now calculate the centrality bins, both 16 & 9
  mcentrality_9 = -1;
  for (int i = 0; i < 9; ++i)
  {
    if (mrefmultcorr >= mcentbin_9[mcentbin_9.size() - i - 1])
    {
      mcentrality_9 = i;
      break;
    }
  }

  mcentrality_16 = -1;
  for (int i = 0; i < 16; ++i) {
    if (mrefmultcorr >= mcentbin_16[mcentbin_16.size() - i - 1]) {
      mcentrality_16 = i;
      break;
    }
  }
  
  // now get the weight
  if(mweight_par.size() && mcentrality_9 >= 0 && mcentrality_16 >= 0 && mrefmultcorr < mweight_bound)
  {
    double par0 = mweight_par[0];
    double par1 = mweight_par[1];
    double par2 = mweight_par[2];
    double par3 = mweight_par[3];
    double par4 = mweight_par[4];
    double par5 = mweight_par[5];
    double par6 = mweight_par[6];
    double ref_const = mrefmultcorr * par2 + par3;
    mweight = par0 + par1 / ref_const + par4 * ref_const + par5 / pow(ref_const, 2.0) + par6 * pow(ref_const, 2.0);
  }
  else
  {
    mweight = 1.0;
  }

  /*
  cout << "refmult " << refmult << ", " << raw_ref << endl;
  cout << "zdc: " << mzdc_par[0] << ", " << mzdc_par[1] << endl;
  cout << "zdc norm: " << mzdc_norm << endl;
  cout << "vz norm: " << mvz_norm << endl;
  cout << "vz: " << mvz_par[0] << ", " << mvz_par[1] << ", " << mvz_par[2] << ", " <<mvz_par[3] << ", " <<mvz_par[4] << ", " <<mvz_par[5] << ", " <<mvz_par[6] << ", " <<endl;
  cout << "ref: " << raw_ref << " -> " << mrefmultcorr << endl;
  cout << "cent bin 16 : " << mcentrality_16 << endl;  
  cout << "weight: " << mweight_par[0] << ", " << mweight_par[1] << ", " << mweight_par[2] << ", " <<mweight_par[3] << ", " <<mweight_par[4] << ", " <<mweight_par[5] << ", " <<mweight_par[6] << ", " <<endl;
  cout << endl;
  */
}
