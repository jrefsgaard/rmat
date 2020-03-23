#ifndef SPECTRUM_DRAWER_H
#define SPECTRUM_DRAWER_H
#include "Spectrum.h"
#include <TH1D.h>

namespace rmat {

class SpectrumDrawer {
  private:
    std::shared_ptr<TH1D> hist;

  public:
    SpectrumDrawer(const char *name, const char *title, int nBins, int binLow, int binHigh);
    ~SpectrumDrawer();
    
    TH1D & Draw(Spectrum &s);
    TH1D & GetHistogram();
};
} //namespace rmat
#endif //SPECTRUM_DRAWER_H
