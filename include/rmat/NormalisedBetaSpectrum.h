#ifndef NORMALISED_BETA_SPECTRUM_H
#define NORMALISED_BETA_SPECTRUM_H
#include <vector>
#include <string>
#include "BetaSpectrum.h"
#include "CompoundSystem.h"

namespace rmat {

class NormalisedBetaSpectrum : public BetaSpectrum {
  /*
  * The NormalisedBetaSpectrum class calculates the spectrum of a beta-delayed 
  * particle decay. Since the norm is known, it is possible to calculate matrix
  * elements and B(GT)'s.
  */
  
  private:
    //The product of N and the halflife can be determined from the total number 
    //of decays or from a specific branch, in which case the partial halflife
    //must be supplied.
    double N;
    double halfLife;
        
  public:
    NormalisedBetaSpectrum(std::string type, CompoundSystem s, double counts = 1., double halflife = 1.);
    ~NormalisedBetaSpectrum();
    
     //Standard order of the parameters are {E,Gamma1, Gamma2,...,B(GT),E,Gamma1,...}.
    virtual void SetParameters(std::vector<double> par);
    virtual void PrintParameters();

    void SetNorm(double counts);
    void SetHalfLife(double halflife);
    void SetMatrixElements(std::vector<double> M);
    void SetMatrixElement(int lambda, double M);
    void SetBGTs(std::vector<double> B);         //BGT = M^2
    void SetBGT(int lambda, double B);

    double GetNorm();
    double GetHalfLife();
    std::vector<double> GetMatrixElements();
    double GetMatrixElement(int lambda);
    std::vector<double> GetBGTs();
    double GetBGT(int lambda);
};
} //namespace rmat;
#endif //NORMALISED_BETA_SPECTRUM_H
