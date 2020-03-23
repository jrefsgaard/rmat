#include "rmat/RMatrixUtils.h"
#include "cpc/cwfcomp.h"
#include <limits>
#include <iostream>
#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_sf_coupling.h>
#include <TMath.h>

using namespace std;
using namespace TMath;

namespace rmat {

double P(int l, double eta, double rho)
{
  /*
  * Returns the penetrability given the angular momentum, l,the Sommerfeld parameter,
  * eta, and rho (k*r).
  */

  /*
  //The implementation with GSL-functions.
  gsl_sf_result F, G;    //The Coulomb functions
  gsl_sf_result Fp, Gp;  //The derivatives.
  double exp_F, exp_G;

  int status = gsl_sf_coulomb_wave_FG_e( eta, rho, l, 0, &F, &Fp, &G, &Gp, &exp_F, &exp_G);
  double fval = F.val;
  double gval = G.val;
  double p = 0.;
  if(status == 16){
    p = rho / (Power(fval,2) * Exp(2 * exp_F) +  Power(gval,2) * Exp(2 * exp_G));
  }
  else{
    p = rho / (Power(fval,2) + Power(gval,2));
  }
  */

  //Implementation with Michel's algorithm. Note that the code has been made unsafe and runs even when result is NaN.
  complex<double> L(l,0);
  complex<double> Eta(eta,0);
  complex<double> z(rho,0);
  class Coulomb_wave_functions cwf(true,L,Eta);
  
  complex<double> F,dF,G,dG;
  cwf.F_dF (z,F,dF);
  cwf.G_dG (z,G,dG);

  //There can be problems with the convergence for very small values of rho. Instead of crashing the calculation we replace
  //infinities or NaNs with largest numbers that can be represented or zero.
  CheckFinity(F,dF,G,dG);

  double p = rho / (Power((double)F.real(),2) + Power((double)G.real(),2));

  return p;
}


double ShiftNeg(int l, double eta, double rho)
{
  /*
  * Returns the shift function given the angular momentum, l, eta (the Sommerfeld parameter)
  * and rho (k*r). For negative energies;
  */

  //The extra factor of 2 is because the differentiation is with respect to rho and not 2*rho.
  double S = rho * 2. *WhittakerWDiff(-eta, l + 0.5, 2 * rho) / WhittakerW(-eta, l + 0.5, 2 * rho);
  return S;
}

double ShiftPos(int l, double eta, double rho)
{
  /*
  * Returns the shift function given the angular momentum, l, eta (the Sommerfeld parameter)
  * and rho (k*r). For positive energies;
  */

  /*
  //Implementation with GSL
  gsl_sf_result F, G;    //The Coulomb functions
  gsl_sf_result Fp, Gp;  //The derivatives.
  double exp_F, exp_G;

  int status = gsl_sf_coulomb_wave_FG_e( eta, rho, l, 0, &F, &Fp, &G, &Gp, &exp_F, &exp_G);
  double S = rho * (F.val * Fp.val + G.val * Gp.val) / (Power(F.val,2) + Power(G.val,2));
  */

  //Implementation with Michel's algorithm.
  complex<double> L(l,0);
  complex<double> Eta(eta,0);
  complex<double> z(rho,0);
  class Coulomb_wave_functions cwf(true,L,Eta);
  
  complex<double> F,dF,G,dG;
  cwf.F_dF (z,F,dF);
  cwf.G_dG (z,G,dG);

  CheckFinity(F,dF,G,dG);

  double S = rho * (F.real() * dF.real() + G.real() * dG.real()) / (Power((double)F.real(),2) + Power((double)G.real(),2));

  return S;
}

double CoulombF(int l, double eta, double rho)
{
  /*
  gsl_sf_result F, G;    //The Coulomb functions
  gsl_sf_result Fp, Gp;  //The derivatives.
  double exp_F, exp_G;

  int status = gsl_sf_coulomb_wave_FG_e( eta, rho, l, 0, &F, &Fp, &G, &Gp, &exp_F, &exp_G);
  double fval = F.val;
  if(status == 16){
    fval = fval * Exp(exp_F);
  }
  */
  complex<double> L(l,0);
  complex<double> Eta(eta,0);
  complex<double> z(rho,0);
  class Coulomb_wave_functions cwf(true,L,Eta);
  
  complex<double> F,dF, G, dG;
  cwf.F_dF (z,F,dF);
  CheckFinity(F,dF,G,dG);

  return F.real();
}

double CoulombG(int l, double eta, double rho)
{
  /*
  gsl_sf_result F, G;    //The Coulomb functions
  gsl_sf_result Fp, Gp;  //The derivatives.
  double exp_F, exp_G;

  int status = gsl_sf_coulomb_wave_FG_e( eta, rho, l, 0, &F, &Fp, &G, &Gp, &exp_F, &exp_G);
  double gval = G.val;
  if(status == 16){
    gval = gval * Exp(exp_G);
  }
  */
  complex<double> L(l,0);
  complex<double> Eta(eta,0);
  complex<double> z(rho,0);
  class Coulomb_wave_functions cwf(true,L,Eta);
  
  complex<double> F, dF, G,dG;
  cwf.G_dG (z,G,dG);
  CheckFinity(F,dF,G,dG);

  return G.real();
}

int CoulombFunctions(int l, double eta, double rho, double &Fval, double &dFval, double &Gval, double &dGval)
{
  /*
  gsl_sf_result F, G;    //The Coulomb functions
  gsl_sf_result Fp, Gp;  //The derivatives.
  double exp_F, exp_G;

  int status = gsl_sf_coulomb_wave_FG_e( eta, rho, l, 0, &F, &Fp, &G, &Gp, &exp_F, &exp_G);

  Fval = F.val;
  Gval = G.val; 

  if(status == 16){
    Fval = Fval * Exp(exp_F);
    Gval = Gval * Exp(exp_G);
  }
  */
  complex<double> L(l,0);
  complex<double> Eta(eta,0);
  complex<double> z(rho,0);
  class Coulomb_wave_functions cwf(true,L,Eta);
  
  complex<double> F,dF,G,dG;
  cwf.F_dF (z,F,dF);
  cwf.G_dG (z,G,dG);
  CheckFinity(F,dF,G,dG);

  Fval = F.real();
  dFval = dF.real();
  Gval = G.real();
  dGval = dG.real();

  return 0;
}

int CoulombFunctions(int l, double eta, double rho, double &Fval, double &Gval)
{
  /*
  gsl_sf_result F, G;    //The Coulomb functions
  gsl_sf_result Fp, Gp;  //The derivatives.
  double exp_F, exp_G;

  int status = gsl_sf_coulomb_wave_FG_e( eta, rho, l, 0, &F, &Fp, &G, &Gp, &exp_F, &exp_G);

  Fval = F.val;
  Gval = G.val; 

  if(status == 16){
    Fval = Fval * Exp(exp_F);
    Gval = Gval * Exp(exp_G);
  }
  */
  double F,dF,G,dG;
  CoulombFunctions(l,eta,rho,F,dF,G,dG);

  Fval = F;
  Gval = G;

  return 0;
}

double WhittakerW(double kappa, double mu, double z)
{
  /*
  * Calculates the Whittaker W-function. The names kappa, mu and z are conventions
  * from Abramowitz & Stegun, and kappa shouldn't be confused with the wave-number.
  */

  double uval = gsl_sf_hyperg_U(mu - kappa + 0.5, 1 + 2*mu, z); //Kummer's confluent hypergeometric U-function.
  double w = Exp(-z/2.) * Power(z, mu + 0.5) * uval;

  return w;
}

//Parameters and GSL-wrapper for the Whittaker-function and its derivative.
struct Wparams {double kappa; double mu;};

double fWhittakerW(double z, void *params)
{
  struct Wparams * p = (struct Wparams *) params;
  return WhittakerW(p->kappa, p->mu, z);
}

double WhittakerWDiff(double kappa, double mu, double z)
{
  /*
  * Calculates the derivative of the Whittaker W-function with respect to z.
  */
  gsl_function W;
  double result, abserr;

  struct Wparams p;
  p.kappa = kappa;
  p.mu = mu;
  W.function = &fWhittakerW;
  W.params = &p;

  gsl_deriv_central (&W, z, 1e-4, &result, &abserr);

  return result;
}

//If the evaluation of the Coulomb functions fail due to the smallness of rho we replace the
//ill defined values with the limiting values for rho=0;
bool CheckFinity(complex<double> &F, complex<double> &dF, complex<double> &G, complex<double> &dG)
{
  bool issue = false;
  //There can be problems with the convergence for very small values of rho. Instead of crashing the calculation we replace
  //infinities or NaNs with largest numbers that can be represented or zero.
  if(!isfinite(F)){F = 0.; issue = true;};
  if(!isfinite(dF)){dF = 0.; issue = true;};
  if(!isfinite(G)){G = numeric_limits<double>::max(); issue = true;};
  if(!isfinite(dG)){dG = -numeric_limits<double>::max(); issue = true;};

  return issue;
}

int Delta(int i, int j)
{
  if(i == j) return 1;
  else return 0;
}

bool IsAlmostZero(double num, double smudgeFactor)
{
  if(abs(num) < smudgeFactor * DBL_EPSILON) return true;
  else return false;
}

double IsAlmostEqual(double num1, double num2)
{
  return IsAlmostZero(num1-num2); 
}

double ClebschGordanFCN(double j1, double j2, double j3, double m1, double m2, double m3)
{
  m3=-m3;
  int j1x2=(int)lrint(2*j1);
  int j2x2=(int)lrint(2*j2);
  int j3x2=(int)lrint(2*j3);
  int m1x2=(int)lrint(2*m1);
  int m2x2=(int)lrint(2*m2);
  int m3x2=(int)lrint(2*m3);

  double w3j=gsl_sf_coupling_3j(j1x2,j2x2,j3x2,m1x2,m2x2,m3x2);

  return pow(-1.0,lrint(j1-j2-m3)) * sqrt(2.0*j3+1.) * w3j;
}

#include <cmath>
#include <complex>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
//#include <gsl/gsl_sf_trig.h>

complex<double> SphericalHarmonic(int l, int m, double theta, double phi)
{
  //cout << "l = " << l << ",  m = " << m << ",  theta = " << theta << ",  phi = " << phi << endl;
  int mabs = abs(m);
  if(mabs > l) return complex<double>(0,0);
  int msign = GSL_SIGN(m);
  double cosT = cos(theta);
  //double cosT = gsl_sf_cos(theta);
  //cout << "cos(theta) = " << cosT << endl;
  double prefix = gsl_sf_legendre_sphPlm(l,mabs,cosT);
  if(msign < 0) prefix *= pow(-1,mabs);// * gsl_sf_fact(l-mabs) / gsl_sf_fact(l+mabs);
  complex<double> i(0,1);
  return prefix * exp(i * m * phi);
  /*
  double sinP, cosP;
  sincosf(m*phi,&sinP,&cosP);
  double re = prefix * cosP;
  double im = prefix * sinP;

  //cout << prefix * exp(i * m * phi) << ",  (" << re << "," << im << ")" << endl;

  return complex<double>(re,im);
  */
}
} //namespace rmat
