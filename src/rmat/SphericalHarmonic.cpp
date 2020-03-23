#include "rmat/SphericalHarmonic.h"
#include <cmath>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <TMath.h>

using namespace std;
namespace rmat {

SphericalHarmonic::SphericalHarmonic(int lmax)
{
  prefactors.resize(lmax+1);
  for(int l=0; l<= lmax; l++){
    prefactors.at(l).resize(2*l+1);
    for(int m=-l; m<=l; m++){
      int mabs = abs(m);
      double factor = sqrt((2.*l+1.)/(4.*TMath::Pi())
                           *gsl_sf_fact(l-mabs)/gsl_sf_fact(l+mabs));
      if(GSL_SIGN(m) < 0){
        factor *= pow(-1,mabs);//*gsl_sf_fact(l-mabs)/gsl_sf_fact(l+mabs);
      }
      prefactors.at(l).at(m+l) = factor;
    }
  }
}

SphericalHarmonic::~SphericalHarmonic(){}

double SphericalHarmonic::P(int l, int m, double x)
{
  //We explicitly implement the first few polynomials.
  if(l == 2){
    if(m == 2) return 3.*(1.-x*x);
    else if(m == 1) return -3.*x*sqrt(1.-x*x);
    else if(m == 0) return 0.5*(3.*x*x - 1.);
    else if(m == -1) return 0.5*x*sqrt(1.-x*x);
    else if(m == -2) return 0.125*(1.-x*x);
  }
  else if(l == 1){
    if(m == 1) return -sqrt(1.-x*x);
    else if(m == 0) return x;
    else if(m == -1) return 0.5*sqrt(1-x*x);
  }
  else if(l == 0) return 1.;
  return gsl_sf_legendre_Plm(l,m,x);
}

complex<double> SphericalHarmonic::Value(int l, int m, double theta, double phi)
{
  int mabs = abs(m);
  if(mabs > l) return complex<double>(0,0);
  double prefactor = prefactors.at(l).at(m+l);
  double cosT = cos(theta);
  double Plm = P(l,mabs,cosT);
  complex<double> i(0,1);
  return prefactor * Plm * exp(i * (double)m * phi);
}
}
