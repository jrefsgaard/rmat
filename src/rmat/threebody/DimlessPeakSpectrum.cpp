#include <rmat/threebody/DimlessPeakSpectrum.h>
#include <armadillo>
#include <iostream>
#include <logft/phase_space.h>

using namespace std;
namespace rmat::threebody {

DimlessPeakSpectrum::DimlessPeakSpectrum(string type, CompoundSystem primSys, CompoundSystem secSys)
: NormalisedBetaSpectrum(type,primSys), secondary_system(secSys)
{
  peakID = 0;
  for(Channel &c : GetSystem().GetChannels()){
    Pair p = c.GetPair();
    int L = c.L();
    double radius = c.Radius();
    Channel chan(p,L);
    chan.SetRadius(radius);
    chan.UseInterpolation();
    channels_2b.push_back(chan);
  }
}

void DimlessPeakSpectrum::SetParameters(vector<double> par)
{
  if(par.size() != NDim()){
    cout << "  DimlessPeakSpectrum::SetParameters(): Wrong number of parameters, par.size() = "<< par.size() << "!" << endl;
    cout << "Expected number: " << NDim() << endl;
    exit(EXIT_FAILURE);
  }
  //First we set the parameters for the primary system.
  int nLevels = GetSystem().GetNLevels();
  int nChannels = GetSystem().GetNChannels();
  //The number of parameters for the primary system.
  int N1 = nLevels * (nChannels + 2);
  for(int l=0; l<nLevels; l++){
    double El = par.at(l*(nChannels + 2));
    double Bl = par.at((l+1)*(nChannels+2)-1);
    std::vector<double> thetal;
    for(int c=0; c<nChannels; c++){
      double thetalc = par.at(l*(nChannels+2) + 1 + c);
      thetal.push_back(thetalc);
    }
    GetSystem().SetEnergy(l,El);
    GetSystem().SetDimensionlessWidths(l,thetal);
    SetBGT(l,Bl);
  }
  
  //Then the secondary compound system.
  int Nlb = secondary_system.GetNLevels();
  int Nc2 = secondary_system.GetNChannels();
  for(int l=0; l<Nlb; l++){
    double El = par.at(N1 + l*(Nc2 + 1));
    vector<double> Gammal;
    for(int c=0; c<Nc2; c++){
      double Gammalc = par.at(N1 + l*(Nc2 + 1) + 1 + c);
      Gammal.push_back(Gammalc);
    }
    //cout << El << "  " << Gammal.at(0) <<  endl;
    secondary_system.SetEnergy(l,El);
    secondary_system.SetPartialWidths(l,Gammal);
  }
}

double DimlessPeakSpectrum::Strength(double Ec)
{
  int A = system.A();
  int Zf = system.Z();
  int Zi = Zf - dZ;
  double W = 0.;
  int Nl = system.GetNLevels(); 
  
  //Calculate renormalisation due to gating on the narrow peak only.
  double renormalisation = secondary_system.GetRenormalisation(peakID);
   
  for(int c : outChannels){
    double Ex = Ec + system.GetThreshold(c);
    //Calculate beta-decay phase space.
    double fBeta = 0.;
    try { fBeta = logft::calculatePhaseSpace(Zi,Zf,A,Ex);} 
    catch (...) { fBeta = 0.;}
    
    //Calculate the penetrability using the standard two-body expression.
    double Pc = channels_2b.at(c).Penetrability(Ec);
    
    vector<double> Js = system.GetJs();
    for(double J : Js){
      //Reduced level matrix for primary system.
      arma::Mat<complex<double>> A = system.LevelMatrix(Ex,J);
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
  
  double strength = W * renormalisation;

  return strength;
}

double DimlessPeakSpectrum::Strength(double Ec, double J)
{
  int A = system.A();
  int Zf = system.Z();
  int Zi = Zf - dZ;
  double W = 0.;
  int Nl = system.GetNLevels();
  
  //Calculate renormalisation due to gating on the narrow peak only.
  double renormalisation = secondary_system.GetRenormalisation(peakID);
  
  for(int c : outChannels){
    double Ex = Ec + system.GetThreshold(c);
    double fBeta = 0.;
    try { fBeta = logft::calculatePhaseSpace(Zi,Zf,A,Ex);} 
    catch (...) { fBeta = 0.;}
    
    double Pc = channels_2b.at(c).Penetrability(Ec);
    
    arma::Mat<complex<double>> A = system.LevelMatrix(Ex,J);
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
  
  double strength = W * renormalisation;
  return strength;
}

void DimlessPeakSpectrum::PrintParameters()
{
  NormalisedBetaSpectrum::PrintParameters();
  cout << "Secondary system:" << endl;
  secondary_system.PrintParameters();
}

int DimlessPeakSpectrum::NDim()
{
  int N2 = secondary_system.GetNLevels() * (secondary_system.GetNChannels() + 1);
  return BetaSpectrum::NDim() + N2;
}

void DimlessPeakSpectrum::SetPeakID(int i)
{
  peakID = i;
}

int DimlessPeakSpectrum::GetPeakID()
{
  return peakID;
}
    
CompoundSystem & DimlessPeakSpectrum::GetSecondarySystem()
{
  return secondary_system;
}
    
vector<Channel> & DimlessPeakSpectrum::GetTwoBodyChannels()
{
  return channels_2b;
}
}
