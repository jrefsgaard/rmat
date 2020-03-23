#include <iostream>
#include <TMath.h>
#include <rmat/RMatrixUtils.h>
#include <rmat/HistSymmetricChiSquare.h>
 
using namespace std;
using namespace TMath;
 
namespace rmat {

HistSymmetricChiSquare::HistSymmetricChiSquare(shared_ptr<TH1> data, shared_ptr<TH1> model) : HistChiSquare{data,model} {}

HistSymmetricChiSquare::~HistSymmetricChiSquare() {}
    
double HistSymmetricChiSquare::Evaluate()
{
  CheckDimension();
  
  double sum = 0.;
  nBins = 0;
  
  int nx = _data->GetXaxis()->GetNbins();
  int ny = _data->GetYaxis()->GetNbins();
  int nz = _data->GetZaxis()->GetNbins();
  
  //Loop over all bins except under-/overflow bins
  for(int zbin=1; zbin<=nz; zbin++){
    for(int ybin=1; ybin<=ny; ybin++){
      for(int xbin=1; xbin<=nx; xbin++){
        double yi = _model->GetBinContent(xbin,ybin,zbin);
        double ni = _data->GetBinContent(xbin,ybin,zbin);
        if(IsAlmostZero(ni + yi)) continue;
        nBins ++;
        sum += Power((ni - yi),2) / (ni + yi);
        //cout << "xbin = " << xbin << ", nbins = " << nBins << ", yi = " << yi << ",  ni = " << ni << endl;
      }
    }
  }
  return 2 * sum;
}
}
