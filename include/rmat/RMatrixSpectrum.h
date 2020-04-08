#ifndef RMATRIX_SPECTRUM_H
#define RMATRIX_SPECTRUM_H
#include "Spectrum.h"
#include "CompoundSystem.h"
#include <vector>

namespace rmat {

class RMatrixSpectrum : public Spectrum {

  protected:
    CompoundSystem system;
    std::vector<int> outChannels;
    
  public:
    RMatrixSpectrum(CompoundSystem s);
    virtual ~RMatrixSpectrum();
    
    virtual double Strength(double Ec) = 0;
    virtual double Strength(double Ec, double J) = 0;
    virtual void SetParameters(std::vector<double> par) = 0;
    virtual void PrintParameters() = 0;
    virtual int NDim() = 0;
    
    void SetOutChannels(std::vector<int> channels);
    void SetOutChannel(int c);

    const std::vector<int> & GetOutChannels();
    CompoundSystem & GetSystem();
};
} //namespace rmat
#endif //RMATRIX_SPECTRUM_H
