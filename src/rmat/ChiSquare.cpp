#include <TMath.h>
#include <iostream>
#include <rmat/ChiSquare.h>
#include <rmat/RMatrixUtils.h>

using namespace std;
using namespace TMath;

namespace rmat {

ChiSquare::ChiSquare(shared_ptr<TH1> d, shared_ptr<Spectrum> m)
{
  data = move(d);
  model = move(m);
  binning = data->GetBinWidth(0);
  minBin = 1;
  maxBin = data->GetNbinsX();
  nBins = 0;
}

ChiSquare::~ChiSquare() {}

double ChiSquare::Evaluate() const
{
  double chi2 = 0.;
  nBins = 0;
  //double counts = data->Integral(minBin,maxBin);
  for(int i=minBin; i<=maxBin; i++){
    double xi = data->GetBinCenter(i);
    //double yi = counts * binning * model->Strength(xi);
    double yi = binning * model->Strength(xi);
    double ni = data->GetBinContent(i);
    double ei = data->GetBinError(i); //Probably sqrt(ni)...?
    if(IsAlmostZero(ei)) ei = 1;
    nBins ++;
    double dChi2 = Power((ni - yi)/ei,2);
    chi2 += dChi2;
  }
  return chi2;
}

void ChiSquare::SetData(shared_ptr<TH1> d)
{
  data = move(d);
}

void ChiSquare::SetModel(shared_ptr<Spectrum> m)
{
  model = move(m);
}

void ChiSquare::SetBinning(double b)
{
  binning = b;
}

void ChiSquare::SetBinRange(int min, int max)
{
  minBin = min;
  maxBin = max;
}

void ChiSquare::SetRange(double min, double max)
{
  minBin = data->FindBin(min);
  maxBin = data->FindBin(max);
}
        
TH1 & ChiSquare::GetData()
{
  return *(data.get());
}

Spectrum & ChiSquare::GetModel()
{
  return *(model.get());
}

double ChiSquare::GetBinning()
{
  return binning;
}

int ChiSquare::GetMinBin()
{
  return minBin;
}

int ChiSquare::GetMaxBin()
{
  return maxBin;
}

double ChiSquare::GetMin()
{
  return data->GetBinLowEdge(minBin);
}

double ChiSquare::GetMax()
{
  return data->GetBinLowEdge(maxBin) + data->GetBinWidth(maxBin);
}

int ChiSquare::GetCounts()
{
  return data->Integral(minBin,maxBin);
}

unsigned int ChiSquare::NDim() const
{
  return model->NDim();
}

int ChiSquare::Ndf()
{
  return GetNBins() - model->NDim();
}

int ChiSquare::GetNBins()
{
  return nBins;
}


double ChiSquare::DoEval(const double *x) const
{
  SetParameters(x);
  return Evaluate();
}

void ChiSquare::SetParameters(const double *par) const
{
  int ndim = model->NDim();
  vector<double> pv;
  pv.assign(par,par + ndim);
  model->SetParameters(pv);
}

ROOT::Math::IBaseFunctionMultiDim * ChiSquare::Clone() const
{
  cout << "ChiSquare::Clone()" << endl;
  ChiSquare *clone = new ChiSquare(*this);  
  return clone;
}
} //namespace rmat
