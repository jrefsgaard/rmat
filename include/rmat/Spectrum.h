#ifndef SPECTRUM_H
#define SPECTRUM_H
//#include <memory>
#include <vector>
#include "CompoundSystem.h"

namespace rmat {

class Spectrum {

  protected:
    //CompoundSystem system;
    //std::vector<int> outChannels;
    //double norm;

  public:
    //Spectrum(CompoundSystem s);
    Spectrum();
    virtual ~Spectrum();
    
    virtual double Strength(double Ec) = 0;
    virtual void SetParameters(std::vector<double> par) = 0;
    virtual void PrintParameters() = 0;
    virtual int NDim() = 0;
    
    //void SetOutChannels(std::vector<int> channels);
    //void SetOutChannel(int c);

    //const std::vector<int> & GetOutChannels();
    //CompoundSystem & GetSystem();
    
    //void SetNorm(double n);
    //double GetNorm();
    
    double Integral(double Emin, double Emax, double epsrel = 1.e-12);
};
} //namespace rmat
#endif //SPECTRUM_H
