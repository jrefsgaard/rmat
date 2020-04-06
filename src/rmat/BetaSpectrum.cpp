#include <iostream>
#include <armadillo>
//#include <stdlib.h>  //exit(), EXIT_FAILURE
#include <TMath.h>
#include <TF1.h>
#include <logft/phase_space.h>
#include <rmat/BetaSpectrum.h>
//#include <rmat/Constants.h>

using namespace std;
using namespace TMath;

namespace rmat {

BetaSpectrum::BetaSpectrum(std::string type, CompoundSystem s) : RMatrixSpectrum(move(s))
{
  SetType(type); 
  Bbeta.resize(system.GetNLevels(),0.);
  //norm = 1.;
}
 
BetaSpectrum::~BetaSpectrum(){}
    
void BetaSpectrum::EnsureSize(std::vector<double> & g)
{
  if(g.size() != system.GetNLevels()){
    cout << "  BetaSpectrum::EnsureSize(): Wrong number of feeding factors. Resizing..." << endl;
    g.resize(system.GetNLevels(),0.);
  }
}
    
double BetaSpectrum::Strength(double Ec)
{
  //cout << "Ec = " << Ec << endl;
  int A = system.A();
  int Zf = system.Z();
  int Zi = Zf - dZ;
  double W = 0.;
  int Nl = system.GetNLevels();  
  for(int c : outChannels){  //Loop over specified outgoing channels.
    double Ex = Ec + system.GetThreshold(c);
    cout << "Ex = " << Ex << endl;
    double fBeta = 0.;
    try { fBeta = logft::calculatePhaseSpace(Zi,Zf,A,Ex);} 
    catch (...) { fBeta = 0.;}                               //Very clever!
    cout << "fBeta = " << fBeta << endl;
    double Pc = system.GetChannel(c).Penetrability(Ec);
    cout << "Pc = " << Pc << endl;
    vector<double> Js = system.GetJs();
    for(double J : Js){
      arma::Mat<complex<double>> A = system.LevelMatrix(Ex,J);  //Reduced level matrix
      //cout << "Level matrix" << endl;
      vector<int> indices = system.GetLevelIndices(J);
      complex<double> amplitude(0,0);
      for(int i=0; i<indices.size(); i++){
        int l = indices.at(i);
        for(int j=0; j<indices.size(); j++){
          int m = indices.at(j);
          double Bl = Bbeta.at(l);
          double gammamc = system.GetWidthAmplitude(m,c);
          amplitude += Bl * gammamc * A(i,j);
        }
      }
      double wc = fBeta * Pc * pow(abs(amplitude),2);  
      W += wc;
    }
  }
  
  double strength = W;
  //double C2 = Log(2)/(Pi()*constants::B);
  //double strength = C2 * W;
  //double strength = halfLife / Log(2.) * C2 * W;
  //double N = 180236.;
  //double halflife = 398.; //Partial halflife, proper halflife is 13.81s
  //double factor =  N * halflife / Log(2.);
  //double binning = 10.;
  //double strength = binning * W * factor * C2; 
  cout << "strength = " << strength << endl;
  return strength;
}

void BetaSpectrum::SetType(std::string type)
{
  if(type == "-") dZ = 1;
  else if(type == "+") dZ = -1;
  else{
    cout << "  BetaSpectrum::SetType(): Unknown decay type: "<< type << "!" << endl;
    exit(EXIT_FAILURE);
  }
}

/*
void BetaSpectrum::SetHalfLife(double t)
{
  if(t < 0.){
    cout << "  BetaSpectrum::SetHalfLife(): Negative half life: "<< t << "!" << endl;
    exit(EXIT_FAILURE);
  }
  halfLife = t;
}
*/

/*
void BetaSpectrum::SetFeedings(std::vector<double> g)
{
  EnsureSize(g);
  gBeta = g;
}

void BetaSpectrum::SetFeeding(int lambda, double g)
{
  if(lambda < 0 || lambda >= gBeta.size()){
    cout << "  BetaSpectrum::SetFeeding(): Index out of range: "<< lambda << "!" << endl;
    exit(EXIT_FAILURE);
  }
  gBeta.at(lambda) = g;
}

void BetaSpectrum::SetMatrixElements(std::vector<double> M)
{
  EnsureSize(M);
  for(int l=0; l<M.size(); l++){
    double Ml = M.at(l);
    SetMatrixElement(l,Ml);
  }
}

void BetaSpectrum::SetMatrixElement(int lambda, double M)
{
  double norm = system.GetRenormalisation(lambda);
  double g = M * Sqrt(norm);
  SetFeeding(lambda,g);
}

void BetaSpectrum::SetBGTs(std::vector<double> B)
{
  EnsureSize(B);
  for(int l=0; l<B.size(); l++){
    double Bl = B.at(l);
    SetBGT(l,Bl);
  } 
}

void BetaSpectrum::SetBGT(int lambda, double B)
{
  //double M = Abs(constants::ga_over_gv) * Sign(Sqrt(Abs(B)),B);
  double M = Sign(Sqrt(Abs(B)),B);
  SetMatrixElement(lambda,M);
}
*/

std::string BetaSpectrum::GetType()
{
  string type;
  if(dZ == 1) type = "-";
  else if(dZ == -1) type = "+";
  return type;
}

/*
double BetaSpectrum::GetHalfLife()
{
  return halfLife;
}
*/

/*
const std::vector<double> BetaSpectrum::GetFeedings()
{
  vector<double> gBeta;
  for(int l=0; l<Bbeta.size(); l++){
    double gl = GetFeeding(l);
    gBeta.push_back(gl);
  }
  return gBeta;
}

double BetaSpectrum::GetFeeding(int lambda)
{
  double Bl = GetB(lambda);
  double gl = Sqrt(Pi() * constants::B / halfLife) * Bl;
  return gl;
}
*/
/*
std::vector<double> BetaSpectrum::GetMatrixElements()
{
  vector<double> M;
  for(int l=0; l<Bbeta.size(); l++){
    double Ml = GetMatrixElement(l);
    M.push_back(Ml);
  }
  return M;
}

double BetaSpectrum::GetMatrixElement(int lambda)
{
  double gl = GetFeeding(lambda);
  double renorm = system.GetRenormalisation(lambda);
  double Ml = gl * Sqrt(renorm);
  return Ml;
}

std::vector<double> BetaSpectrum::GetBGTs()
{
  vector<double> B;
  for(int l=0; l<Bbeta.size(); l++){
    double Bl = GetBGT(l);
    B.push_back(Bl);
  }
  return B;
}

double BetaSpectrum::GetBGT(int lambda)
{
  double Ml = GetMatrixElement(lambda);
  //double Bl = Power(constants::ga_over_gv,-2) * Power(Ml,2);
  double Bl = Sign(Power(Ml,2),Ml);
  return Bl;
}
*/

void BetaSpectrum::PrintParameters()
{
  cout << "--Beta-decay parameters--" << endl;
  //cout << "  Calculated for N = " << norm << endl;
  int N = Bbeta.size();
  for(int i=0; i<N; i++){
    cout << "  Level " << i << ":" << endl;
    //cout << "    g = " << GetFeeding(i) << endl;
    cout << "    B = " << GetB(i) << endl;
    //cout << "    M = " << GetMatrixElement(i) << endl;
    //cout << "    B(GT) = " << GetBGT(i) << endl;
  }
  cout << endl;
}

//How many parameters does the spectrum depend on?
int BetaSpectrum::NDim()
{
  //The spectrum depends on N = Nlevels * (Nchannels + 2)
  //The +2 is for the level energy and the beta decay feeding factor.
  int nLevels = GetSystem().GetNLevels();
  int nChannels = GetSystem().GetNChannels();
  return nLevels * (nChannels + 2);
}

void BetaSpectrum::SetParameters(vector<double> par)
{
  if(par.size() != NDim()){
    cout << "  BetaSpectrum::SetParameters(): Wrong number of parameters, par.size() = "<< par.size() << "!" << endl;
    exit(EXIT_FAILURE);
  }
  int nLevels = GetSystem().GetNLevels();
  int nChannels = GetSystem().GetNChannels();
  for(int l=0; l<nLevels; l++){
    double El = par.at(l*(nChannels + 2));
    double Bl = par.at((l+1)*(nChannels+2)-1);
    vector<double> Gammal;
    for(int c=0; c<nChannels; c++){
      double Gammalc = par.at(l*(nChannels+2) + 1 + c);
      Gammal.push_back(Gammalc);
    }
    //cout << El << "  " << Gammal.at(0) << "  " << Bl << endl;
    GetSystem().SetEnergy(l,El);
    GetSystem().SetPartialWidths(l,Gammal);
    SetB(l,Bl);
  }
}

void BetaSpectrum::SetBs(std::vector<double> B)
{
  EnsureSize(B);
  Bbeta = B;
}

void BetaSpectrum::SetB(int lambda, double B)
{
  if(lambda < 0 || lambda >= Bbeta.size()){
    cout << "  BetaSpectrum::SetB(): Index out of range: "<< lambda << "!" << endl;
    exit(EXIT_FAILURE);
  }
  Bbeta.at(lambda) = B;
}

std::vector<double> & BetaSpectrum::GetBs()
{
  return Bbeta;  
}

double BetaSpectrum::GetB(int lambda)
{
  return Bbeta.at(lambda);
}

/*
void BetaSpectrum::SetNorm(double N)
{
  if(N < 0.){
    cout << "  BetaSpectrum::SetNorm(): Invalid value of argument: "<< N << "!" << endl;
    exit(EXIT_FAILURE);
  }
  norm = N;
}

double BetaSpectrum::GetNorm()
{
  return norm;
}
*/
} //namespace rmat;
