#include <rmat/Spectrum.h>
#include <TF1.h>

using namespace std;

namespace rmat {

//Spectrum::Spectrum(CompoundSystem s) : system(move(s)) {norm = 1.;}
//Spectrum::Spectrum() {norm = 1.;}
Spectrum::Spectrum() {}

Spectrum::~Spectrum(){}

/*
void Spectrum::SetOutChannels(std::vector<int> channels)
{
  int Nc = system.GetNChannels();
  for(int ci : channels){
    if(ci < 0 || ci >= Nc){
      cout << "  Spectrum::SetOutChannels(): Channel " << ci << " does not exist!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  outChannels = channels;
}

void Spectrum::SetOutChannel(int c)
{
  vector<int> cv = {c};
  SetOutChannels(cv);
}

const vector<int> & Spectrum::GetOutChannels()
{
  return outChannels;
}

CompoundSystem & Spectrum::GetSystem()
{
  return system;
}


void Spectrum::SetNorm(double N)
{
  if(N < 0.){
    cout << "  Spectrum::SetNorm(): Invalid value of argument: "<< N << "!" << endl;
    exit(EXIT_FAILURE);
  }
  norm = N;
}

double Spectrum::GetNorm()
{
  return norm;
}
*/

double Spectrum::Integral(double Emin, double Emax, double epsrel)
{
  auto spectrum = [&](const double* x, const double *p) mutable {
    return Strength(x[0]);
  };
  TF1 f("fIntegrand",[&](double *X, double *P){ return spectrum(X,P);},Emin,Emax,0);
    
  return f.Integral(Emin,Emax,epsrel);
}
} //namespace rmat;
