#include <stdlib.h>
#include <iostream>
#include <vector>
#include <complex>
#include <armadillo>
#include <TMath.h>
#include <rmat/threebody/DensityOfStates.h>

using namespace std;

namespace rmat::threebody {

DensityOfStates::DensityOfStates(int A, int Z) : system {A,Z,1} {}

DensityOfStates::~DensityOfStates() {}

void DensityOfStates::SetLevel(double E, double Gamma)
{
  if(system.GetNChannels() == 0){
    cout << "  DensityOfStates::SetLevel(): Channels must be specified before the partial width can be set!" << endl;
    exit(EXIT_FAILURE);
  }
  system.ClearLevels();
  system.AddLevel(E);
  vector<double> partialWidth {Gamma};
  system.SetPartialWidths(0,partialWidth);
}

void DensityOfStates::SetChannel(Channel c, double threshold)
{
  system.ClearChannels();
  system.AddChannel(move(c),threshold);
}
 
void DensityOfStates::SetChannel(std::string type, int l, double r0, double threshold)
{
  system.ClearChannels();
  system.AddChannel(type,l,r0,threshold);
}

void DensityOfStates::UseInterpolation(bool interp, double Emax, double Estep)
{
  system.GetChannel(0).UseInterpolation(interp,Emax,Estep);
}

double DensityOfStates::Evaluate(double Ec)
{
  double Ex = Ec + system.GetThreshold(0);
  vector<double> Js = system.GetJs(); //Here, irrelevant but necessary.
  double J = Js.at(0);
  arma::Mat<complex<double>> A = system.LevelMatrix(Ex,J);
  complex<double> amplitude(0,0);
  double gammamc = system.GetWidthAmplitude(0,0);
  amplitude += gammamc * A(0,0);
  double Pc = system.GetChannel(0).Penetrability(Ec);
  double rho = Pc * pow(abs(amplitude),2);
  
  return rho / TMath::Pi();
}
}//namespace rmat::threebody
