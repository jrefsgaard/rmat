#include <limits>
#include <iostream>
#include <stdlib.h>  //exit(), EXIT_FAILURE
#include <armadillo>
#include <TMath.h>
#include <rmat/CompoundSystem.h>
#include <rmat/Nucleus.h>
#include <rmat/Pair.h>
#include <rmat/RMatrixUtils.h>

using namespace std;
using namespace arma;
using namespace TMath;

namespace rmat {

CompoundSystem::CompoundSystem(int A, int Z, int Nl) : a(A), z(Z)
{
  for(int i=0; i<Nl; i++){
    Level li(0.,vector<double>());
    levels.push_back(li);
  }
}

CompoundSystem::~CompoundSystem(){}

void CompoundSystem::CheckAndResize()
{
  for(Level & li : levels){
    EnsureWidthSize(li.GetWidthAmplitudes());
  }
}

void CompoundSystem::EnsureWidthSize(vector<double> & widths)
{
  int Nc = GetNChannels();
  widths.resize(Nc,0);
}

Channel & CompoundSystem::AddChannel(Channel c, double threshold)
{
  channels.push_back(move(c));
  thresholds.push_back(threshold);
  CheckAndResize();
  return channels.back();
}

Channel & CompoundSystem::AddChannel(string type, int l, double r0, double threshold)
{
  if(type == "a"){
    //The user wants to add an alpha decay channel.
    int Ares = A() - 4;
    int Zres = Z() - 2;
    Nucleus alpha_particle(4,2);
    Nucleus residual(Ares,Zres);
    Pair p(residual,alpha_particle);
    Channel c(p,l,r0);
    return AddChannel(move(c),threshold);
  }
  else if(type == "p"){
    //Proton channel.
    int Ares = A() - 1;
    int Zres = Z() - 1;
    Nucleus proton(1,1);
    Nucleus residual(Ares,Zres);
    Pair p(residual,proton);
    Channel c(p,l,r0);
    return AddChannel(move(c),threshold);
  }
  else if(type == "d"){
    //Deuteron channel.
    int Ares = A() - 2;
    int Zres = Z() - 1;
    Nucleus deuteron(2,1);
    Nucleus residual(Ares,Zres);
    Pair p(residual,deuteron);
    Channel c(p,l,r0);
    return AddChannel(move(c),threshold);
  }
  else{
    cout << "  CompoundSystem::AddChannel(): Unknown channel type: "<< type << "!" << endl;
    exit(EXIT_FAILURE);
  }
}

void CompoundSystem::AddLevel(Level level)
{
  levels.push_back(level);
  CheckAndResize();
}

void CompoundSystem::AddLevel(double E, std::vector<double> gamma)
{
  Level l(E, gamma);
  AddLevel(l);
}
 
void CompoundSystem::AddLevel(double E)
{
  Level l(E, vector<double>());
  AddLevel(l);
}

void CompoundSystem::ClearLevels(){ levels.clear();}

void CompoundSystem::ClearChannels()
{
  channels.clear();
  thresholds.clear();
  CheckAndResize();
}

void CompoundSystem::UseInterpolation(bool interp, double Emax, double Estep)
{
  int Nc = GetNChannels();
  for(int c=0; c<Nc; c++){
    GetChannel(c).UseInterpolation(interp,Emax,Estep);
  }
}

void CompoundSystem::SetEnergy(int lambda, double E)
{
  //cout << "SetEnergy(): l = " << lambda << ", E = " << E << endl;
  GetLevel(lambda).SetE(E);
}

void CompoundSystem::SetPartialWidths(int lambda, std::vector<double> Gamma)
{
  EnsureWidthSize(Gamma);
  double El = GetEnergy(lambda);
  vector<double> gamma2ObsVector;
  double sum = 0.;
  for(int c=0; c<Gamma.size(); c++){
    Channel & channel = GetChannel(c);
    double Gammac = Gamma.at(c);
    double Ec = GetChannelEnergy(El,c);
    double Pc = channel.Penetrability(Ec);
    double dSc = channel.ShiftFunctionDeriv(Ec);
    //cout << Ec << "  " << Pc << "  " << dSc << "  ";
    double gamma2Obs = Abs(Gammac / (2. * Pc));
    //cout << gamma2Obs << endl;
    gamma2ObsVector.push_back(gamma2Obs);
    sum += gamma2Obs * dSc;
  }
  if(sum > 1.){
    cout << "CompoundSystem: Warning! Requested width is too large. Expect unphysical output." << endl;
  }
  for(int c=0; c<Gamma.size(); c++){
    double gamma2 = gamma2ObsVector.at(c) / (1. - sum);
    //cout << "  " << gamma2 << endl;
    SetReducedWidth(lambda,c,Sign(gamma2,Gamma.at(c)));
  }
}

void CompoundSystem::SetWidthAmplitude(int lambda, int c, double gamma)
{
  GetLevel(lambda).SetWidthAmplitude(c,gamma);
}

void CompoundSystem::SetWidthAmplitudes(int lambda, std::vector<double> gamma)
{
  EnsureWidthSize(gamma);
  GetLevel(lambda).SetGammas(gamma);
}

void CompoundSystem::SetReducedWidth(int lambda, int c, double gamma2)
{
  double gamma = Sign(Sqrt(Abs(gamma2)),gamma2);
  SetWidthAmplitude(lambda,c,gamma);
}

void CompoundSystem::SetReducedWidths(int lambda, std::vector<double> gamma2)
{
  vector<double> gamma;
  for(int c=0; c<gamma2.size(); c++){
    double g = Sign(Sqrt(Abs(gamma2.at(c))),gamma2.at(c));
    gamma.push_back(g);
  }
  SetWidthAmplitudes(lambda,gamma);
}

void CompoundSystem::SetDimensionlessWidth(int lambda, int c, double theta2)
{
  double limit = GetChannel(c).WignerLimit();
  double gamma2 = 2./3. * limit * theta2;
  SetReducedWidth(lambda,c,gamma2);
}

void CompoundSystem::SetDimensionlessWidths(int lambda, vector<double> theta2)
{
  vector<double> gamma2;
  for(int c=0; c<theta2.size(); c++){
    double limit = GetChannel(c).WignerLimit();
    double g2 = 2./3. * limit * theta2.at(c);
    gamma2.push_back(g2);
  }
  SetReducedWidths(lambda,gamma2);
}

double CompoundSystem::GetEnergy(int lambda)
{
  double El = GetLevel(lambda).E();
  return El;
}

double CompoundSystem::GetChannelEnergy(double Ex, int c)
{
  double Ec = Ex - GetThreshold(c);
  return Ec;
}

double CompoundSystem::GetPartialWidth(int lambda, int c)
{
  double gamma2 = GetReducedWidth(lambda,c);
  double El = GetEnergy(lambda);
  double Ec = GetChannelEnergy(El,c);
  double Pc = GetChannel(c).Penetrability(Ec);
  double renormalisation = GetRenormalisation(lambda);
  double Gamma = 2. * Pc * gamma2 * renormalisation;
  return Gamma;
}

double CompoundSystem::GetWidthAmplitude(int lambda, int c)
{
  return GetLevel(lambda).GetWidthAmplitude(c);
}

double CompoundSystem::GetReducedWidth(int lambda, int c)
{
  double gamma = GetWidthAmplitude(lambda,c);
  return Power(gamma,2);
}

double CompoundSystem::GetDimensionlessWidth(int lambda, int c)
{
  double gammalc = GetReducedWidth(lambda,c);
  double limit = GetChannel(c).WignerLimit();
  return gammalc / (2./3. * limit);
}

double CompoundSystem::GetWidth(int lambda)
{
  double width = 0.;
  for(int c=0; c<GetNChannels(); c++){
    width += GetPartialWidth(lambda,c);
  }
  return width;
}

double CompoundSystem::GetRenormalisation(int lambda)
{
  double sum = 1.;
  int Nc = GetNChannels();
  double El = GetEnergy(lambda);
  for(int ci=0; ci<Nc; ci++){
    double Ec = GetChannelEnergy(El,ci);
    sum += GetReducedWidth(lambda,ci) * GetChannel(ci).ShiftFunctionDeriv(Ec);
  }
  return 1. / sum;
}

Level & CompoundSystem::GetLevel(int l)
{
  return levels.at(l);
}

vector<Level> & CompoundSystem::GetLevels()
{
  return levels;
}

Channel & CompoundSystem::GetChannel(int c)
{
  return channels.at(c);
}

vector<Channel> & CompoundSystem::GetChannels()
{
  return channels;
}

double CompoundSystem::GetThreshold(int c)
{
  return thresholds.at(c);
}

vector<double> & CompoundSystem::GetThresholds()
{
  return thresholds;
}

int CompoundSystem::GetNLevels(){ return levels.size();}

int CompoundSystem::GetNChannels(){ return channels.size();}

int CompoundSystem::A(){return a;}

int CompoundSystem::Z(){return z;}

//Returns diminished level matrix
Mat<complex<double>> CompoundSystem::LevelMatrix(double Ex, double J)
{
  //cout << "line 273 J = " << J << endl;
  vector<int> indices = GetLevelIndices(J);
  Mat<complex<double>> A(indices.size(),indices.size());
  
  //if(indices.size() != 0) cout << "indices " << indices.at(0) << endl;
  
  complex<double> i(0,1);
  for(int j=0; j<indices.size(); j++){
    //cout << "line 279" << endl;
    int l = indices.at(j);
    double El = GetEnergy(l);
    //cout << "line 282" << endl;
    for(int k=0; k<=j; k++){
      //cout << "line 284" << endl;
      int m = indices.at(k);
      //cout << "line 286" << endl;
      double Em = GetEnergy(m);
      //cout << "line 288" << endl;
      complex<double> value = (El - Ex) * Delta(l,m);
      for(int c=0; c<channels.size(); c++){
        //cout << "line 291" << endl;
        Channel & channel = GetChannel(c);
        double Ethres = GetThreshold(c);
        double Ec = Ex - Ethres;
        //cout << "Ex = " << Ex << ", Ethres = " << Ethres << ", Ec = " << Ec << endl;
        double gammalc = GetWidthAmplitude(l,c);
        double gammamc = GetWidthAmplitude(m,c);
        double Sc = channel.ShiftFunction(Ec);
        double Pc = channel.Penetrability(Ec);
        value -= gammalc * gammamc * (Sc + i * Pc);
        //cout << "line 33" << endl;
        if(l == m){
          //cout << "El - Ethresh = " << El - Ethres << endl;
          double Slc = channel.ShiftFunction(El - Ethres);
          //cout << "line 336" << endl;
          value += pow(gammalc,2) * Slc;
        }
        else{
          double Slc = channel.ShiftFunction(El - Ethres);
          double Smc = channel.ShiftFunction(Em - Ethres);
          value += gammalc * gammamc * (Slc * (Ex - Em) - Smc * (Ex - El)) / (El - Em);
        }
      }
      A(j,k) = value;
    }
  }
  A = symmatl(A,false); //Reflect lower triangle into the upper triangle.
    
  return inv(A);
}

const void CompoundSystem::PrintParameters()
{
  int Nl = GetNLevels();
  int Nc = GetNChannels();
  cout << "--On-resonance R-matrix parameters--" << endl ;
  for(int l=0; l<Nl; l++){
    cout << "  Level " << l << ":" << endl;
    cout << "    E = " << GetEnergy(l) << endl;
    cout << "    Gamma = " << GetWidth(l) << endl;
    for(int c=0; c<Nc; c++){
      cout << "    g(" << l << c << ") = " << GetWidthAmplitude(l,c) ;
      cout << ",  Gamma(" << l << c << ") = " << GetPartialWidth(l,c) ;
      cout << ",  theta(" << l << c << ") = " << GetDimensionlessWidth(l,c) << endl;
    }
  }
  cout << endl ;
}

vector<double> CompoundSystem::GetJs()
{
  vector<double> Js;
  for(Level l : levels){
    double Jl = l.J();
    bool hasJ = false;
    for(double J : Js){
      if(IsAlmostEqual(J,Jl)) hasJ = true;
    }
    if(!hasJ) Js.push_back(Jl);
  } 
  return Js;
}

vector<int> CompoundSystem::GetLevelIndices(double J)
{
  //cout << "line 361" << endl;
  vector<int> indices;
  for(int l=0; l<GetNLevels(); l++){
    //cout << "l = " << l << ",  GetNLevels() = " << GetNLevels() << endl;
    double Jl = levels.at(l).J();
    if(IsAlmostEqual(Jl,J)){
      indices.push_back(l);
      //cout << "level " << l << " has the right spin, Jl = " << Jl << " compared to J = " << J << endl;
    }
  }
  return indices;
}

} //namespace rmat
