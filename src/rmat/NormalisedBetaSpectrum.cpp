#include <iostream>
#include <TMath.h>
#include <rmat/NormalisedBetaSpectrum.h>
#include <rmat/Constants.h>

using namespace std;
using namespace TMath;

namespace rmat {

NormalisedBetaSpectrum::NormalisedBetaSpectrum(std::string type, CompoundSystem s, double counts, double halflife)
: BetaSpectrum(type, move(s)) , N(counts), halfLife(halflife)
{
  N = counts;
  halfLife = halflife;
}

NormalisedBetaSpectrum::~NormalisedBetaSpectrum() {}
    
void NormalisedBetaSpectrum::SetNorm(double counts)
{
  if(counts < 0.){
    cout << "  BetaSpectrum::SetNorm(): Invalid value of argument: "<< counts << "!" << endl;
    exit(EXIT_FAILURE);
  }
  N = counts;
}

void NormalisedBetaSpectrum::SetHalfLife(double halflife)
{
  if(halflife < 0.){
    cout << "  BetaSpectrum::SetHalfLife(): Negative half life: "<< halflife << "!" << endl;
    exit(EXIT_FAILURE);
  }
  halfLife = halflife;
}

void NormalisedBetaSpectrum::SetMatrixElements(vector<double> M)
{
  EnsureSize(M);
  for(int l=0; l<M.size(); l++){
    double Ml = M.at(l);
    SetMatrixElement(l,Ml);
  }
}

void NormalisedBetaSpectrum::SetMatrixElement(int lambda, double M)
{
  double norm = system.GetRenormalisation(lambda); // = 1/sqrt(1 + ...)
  double B = M * Sqrt(N * halfLife /(Pi() * constants::B)) * Sqrt(1./norm);
  SetB(lambda,B);
}

void NormalisedBetaSpectrum::SetBGTs(vector<double> B) //BGT = M^2
{
  EnsureSize(B);
  for(int l=0; l<B.size(); l++){
    double Bl = B.at(l);
    SetBGT(l,Bl);
  } 
}

void NormalisedBetaSpectrum::SetBGT(int lambda, double B)
{
  double M = Abs(constants::ga_over_gv) * Sign(Sqrt(Abs(B)),B);
  //double M = Sign(Sqrt(Abs(B)),B);
  SetMatrixElement(lambda,M);
}

double NormalisedBetaSpectrum::GetNorm()
{
  return N;
}

double NormalisedBetaSpectrum::GetHalfLife()
{
  return halfLife;
}

vector<double> NormalisedBetaSpectrum::GetMatrixElements()
{
  vector<double> M;
  for(int l=0; l<Bbeta.size(); l++){
    double Ml = GetMatrixElement(l);
    M.push_back(Ml);
  }
  return M;
}

double NormalisedBetaSpectrum::GetMatrixElement(int lambda)
{
  double Bl = GetB(lambda);
  double renorm = system.GetRenormalisation(lambda);
  double Ml = Sqrt(Pi() * constants::B / (N * halfLife)) * Sqrt(renorm) * Bl;
  return Ml;
}

vector<double> NormalisedBetaSpectrum::GetBGTs()
{
  vector<double> B;
  for(int l=0; l<Bbeta.size(); l++){
    double Bl = GetBGT(l);
    B.push_back(Bl);
  }
  return B;
}

double NormalisedBetaSpectrum::GetBGT(int lambda)
{
  double Ml = GetMatrixElement(lambda);
  double Bl = Sign(Power(Ml,2),Ml) * Power(1./(constants::ga_over_gv),2);
  return Bl;
}
     //Standard order of the parameters are {E,Gamma1, Gamma2,...,B(GT),E,Gamma1,...}.
void NormalisedBetaSpectrum::SetParameters(vector<double> par)
{
  BetaSpectrum::SetParameters(par);
  //Actually, the user meant to supply B(GT)'s...
  vector<double> BGT = GetBs();
  SetBGTs(BGT);
}

void NormalisedBetaSpectrum::PrintParameters()
{
  cout << "--Beta-decay parameters--" << endl;
  //cout << "  Calculated for N = " << norm << endl;
  int N = Bbeta.size();
  for(int i=0; i<N; i++){
    cout << "  Level " << i << ":" << endl;
    cout << "    B = " << GetB(i) << endl;
    cout << "    M = " << GetMatrixElement(i) << endl;
    cout << "    B(GT) = " << GetBGT(i) << endl;
  }
  cout << endl;
}
} //namespace rmat;
