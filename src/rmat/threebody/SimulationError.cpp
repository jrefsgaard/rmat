#include <iostream>
#include <TH3.h>
#include <TH2.h>
#include "rmat/threebody/SimulationError.h"
//#include "rmat/threebody/Simulator.h"
#include "rmat/HistLikelihoodChiSquare.h"
#include <TFile.h>

using namespace std;
namespace rmat::threebody {

SimulationError::SimulationError(Simulator model, std::shared_ptr<TH1> data)
: _model(model), _data(data) {normalise = true;}

SimulationError::~SimulationError(){}

double SimulationError::DoEval(const double *x) const
{
  SetParameters(x);
  return Evaluate();
}

double SimulationError::Evaluate() const
{
  //Perform simulation.
  long long int N = _model.Simulate();
  //cout << "N = " << N << endl;
  
  //Obtain simulated spectrum.
  int nx = _data->GetXaxis()->GetNbins();
  double xmin = _data->GetXaxis()->GetXmin();
  double xmax = _data->GetXaxis()->GetXmax();
  //cout << "nx = " << nx << ", xmin = " << xmin << ", xmax = " << xmax << endl;
  int ny = _data->GetYaxis()->GetNbins();
  double ymin = _data->GetYaxis()->GetXmin();
  double ymax = _data->GetYaxis()->GetXmax();
  //cout << "ny = " << ny << ", ymin = " << ymin << ", ymax = " << ymax << endl;
  int nz = _data->GetZaxis()->GetNbins();
  double zmin = _data->GetZaxis()->GetXmin();
  double zmax = _data->GetZaxis()->GetXmax();
  //cout << "nz = " << nz << ", zmin = " << zmin << ", zmax = " << zmax << endl;
  shared_ptr<TH1> result = _model.GetSpectrum(nx,xmin,xmax,ny,ymin,ymax,nz,zmin,zmax);

  //TFile outfile("~/workspace/i161/fit_test/test.root","RECREATE");
  //_data->Write("data");
  //result->Write("result");
  
  if(normalise){
    double integral;
    //We calculate the integral of the data spectrum, remembering excluded Q-values.
    if(nx > 1 && ny > 1 && nz > 1){
      //Standard 3D-spectrum with (x,y,Q)-coordinates.
      auto cast_histogram = dynamic_pointer_cast<TH3>(_data);
      integral = cast_histogram->Integral();
      //cout << "integral = " << integral << endl;
      vector<array<double,2>> ranges = _model.GetDecayWeight().GetExcludedQRanges();
      for(auto range : ranges){
        //These limits assume that the user ranges correspond to whole bins.
        int minBin = cast_histogram->GetZaxis()->FindBin(range[0]);
        int maxBin = cast_histogram->GetZaxis()->FindBin(range[1])-1;
        integral -= cast_histogram->Integral(1,nx,1,ny,minBin,maxBin);
        //cout << "minBin = " << minBin << ",  maxBin = " << maxBin << ", integral = " << cast_histogram->Integral(1,nx,1,ny,minBin,maxBin) << endl;
      }
    }
    else if(nx > 1 && ny > 1 && nz <= 1){
      //Standard 2D-spectrum with (x,y)-coordinates.
      auto cast_histogram = dynamic_pointer_cast<TH2>(_data);
      integral = cast_histogram->Integral();
    }
    else if(nx > 1 && ny <= 1 && nz <= 1){
      //Standard 1D-spectrum with Q-values.
      integral = _data->Integral();
      //cout << "integral = " << integral << endl;
      vector<array<double,2>> ranges = _model.GetDecayWeight().GetExcludedQRanges();
      for(auto range : ranges){
        int minBin = _data->GetXaxis()->FindBin(range[0]);
        int maxBin = _data->GetXaxis()->FindBin(range[1])-1;
        integral -= _data->Integral(minBin,maxBin);
        //cout << "minBin = " << minBin << ",  maxBin = " << maxBin << ", integral = " << _data->Integral(minBin,maxBin) << endl;
      }     
    }
    else{
      cout << "  SimulationError::Evaluate(): Unknown data binning:" ;
      cout << " nx = " << nx << ", ny = " << ny << ", nz = " << nz << endl;
      exit(EXIT_FAILURE);
    }
    double scale = integral / N;
    result->Scale(scale);
  }

  //result->Write("scaled_result");
  //outfile.Close();  

  //Calculate error.
  HistLikelihoodChiSquare error(_data,result);
  double errValue = error.Evaluate();

  nBins = error.GetNBins();

  //cout << "Error = " << errValue << ",  nBins = " << nBins << ",  Nsim = " << N << ",  Ndata = " << _data->Integral() << endl;
  
  return errValue;
}

unsigned int SimulationError::NDim() const
{
  return _model.NDim();
}

void SimulationError::SetParameters(const double *par) const
{
  int ndim = _model.NDim();
  vector<double> pv;
  pv.assign(par,par + ndim);
  _model.SetParameters(pv);
  //XXX
  //_model.GetDecayWeight().PrintParameters();
}

void SimulationError::PrintParameters()
{
  _model.GetDecayWeight().PrintParameters();
}

ROOT::Math::IBaseFunctionMultiDim * SimulationError::Clone() const
{
  SimulationError *clone = new SimulationError(*this);
  return clone;
}

void SimulationError::SetModel(Simulator model)
{
  _model = model;
}

void SimulationError::SetData(std::shared_ptr<TH1> data)
{
  _data = data;
}

Simulator & SimulationError::GetModel()
{
  return _model;
}

shared_ptr<TH1> SimulationError::GetData()
{
  return _data;
}

int SimulationError::Ndf()
{
  return GetNBins() - NDim();
}

int SimulationError::GetNBins() const
{
  //cout << "SimulationError::GetNBins(): nBins = " << nBins << endl;
  return nBins;
}

void SimulationError::AutoNormalise(bool arg)
{
  normalise = arg;
}
}//namespace rmat::threebody
