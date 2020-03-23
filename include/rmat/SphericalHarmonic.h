#ifndef SPHERICAL_HARMONIC_H
#define SPHERICAL_HARMONIC_H
#include <complex>
#include <vector>

namespace rmat {

class SphericalHarmonic {
  private:
    std::vector<std::vector<double>> prefactors;

    /**
    * Explicitly implements the first three orders of the Legendre
    * polynomials. Higher orders are calculated with GSL.
    */
    double P(int l, int m, double x);
    
  public:
    SphericalHarmonic(int lmax = 10);
    ~SphericalHarmonic();
    
    std::complex<double> Value(int l, int m, double theta, double phi);
};
}
#endif
