#include <iostream>
#include <stdlib.h>
#include <TMath.h>
#include <rmat/RMatrixUtils.h>
#include <rmat/HistChiSquare.h>

using namespace std;
using namespace TMath;

namespace rmat {
    
HistChiSquare::HistChiSquare(shared_ptr<TH1> data, shared_ptr<TH1> model) : _data(data), _model(model) {}

HistChiSquare::~HistChiSquare() {}
    
void HistChiSquare::SetData(shared_ptr<TH1> data)
{
  _data = data;
}

void HistChiSquare::SetModel(shared_ptr<TH1> model)
{
  _model = model;
}

TH1 & HistChiSquare::GetData()
{
  return *(_data.get());
}

TH1 & HistChiSquare::GetModel()
{
  return *(_model.get());
}
 
int HistChiSquare::GetNBins()
{
  //cout << "HistChiSquare::GetNBins(): nBins = " << nBins << endl;
  return nBins;
}

void HistChiSquare::CheckDimension()
{
  if(_data->GetXaxis()->GetNbins() != _model->GetXaxis()->GetNbins()){
    cout << "  HistChiSquare::CheckDimension(): Data and model histogram have different x-dimension!" << endl;
    exit(EXIT_FAILURE);
  }
  if(_data->GetYaxis()->GetNbins() != _model->GetYaxis()->GetNbins()){
    cout << "  HistChiSquare::CheckDimension(): Data and model histogram have different y-dimension!" << endl;
    exit(EXIT_FAILURE);
  }
  if(_data->GetZaxis()->GetNbins() != _model->GetZaxis()->GetNbins()){
    cout << "  HistChiSquare::CheckDimension(): Data and model histogram have different z-dimension!" << endl;
    exit(EXIT_FAILURE);
  }
}
   
double HistChiSquare::Evaluate()
{
  CheckDimension();
  
  double chi2 = 0.;
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
        double ei = _data->GetBinError(xbin,ybin,zbin); //Probably sqrt(ni)...?
        if(IsAlmostZero(ei)) ei = 1;
        nBins ++;
        double dChi2 = Power((ni - yi)/ei,2);
        chi2 += dChi2;
      }
    }
  }
  return chi2;
}
}
