#ifndef BETA_SPECTRUM_H
#define BETA_SPECTRUM_H
#include <vector>
#include <memory>
#include <string>
#include "Spectrum.h"
#include "RMatrixSpectrum.h"
#include "CompoundSystem.h"

namespace rmat {

class BetaSpectrum : public RMatrixSpectrum {
  /*
  * The BetaSpectrum class calculates the spectrum of a beta-delayed particle
  * decay.
  */
  
  private:
    int dZ;   //1 (beta-minus) or -1 (beta-plus)

  protected: 
    std::vector<double> Bbeta;
   
    void EnsureSize(std::vector<double> & g);  //Ensures that the vector has proper length.
    
  public:
    BetaSpectrum(std::string type, CompoundSystem s);
    ~BetaSpectrum();
    
    //These methods MUST be implemented.
    double Strength(double Ec);
    //Standard order of the parameters are {E,Gamma1, Gamma2,...,B,E,Gamma1,...}.
    virtual void SetParameters(std::vector<double> par);
    int NDim();
    virtual void PrintParameters();
    
    void SetType(std::string type);
    //There are at least four different ways to characterise the
    //beta strength to a given level.
    //We use the B(lambda)'s from Barker and Warburton
    void SetBs(std::vector<double> B);
    void SetB(int lambda, double B);
    
    std::string GetType();
    std::vector<double> & GetBs();
    double GetB(int lambda);
};
} //namespace rmat;
#endif //BETA_SPECTRUM_H
