#include "rmat/ClebschGordan.h"
#include <cmath>
#include <assert.h>
#include <iostream>
#include <gsl/gsl_sf_coupling.h>
#include <TMath.h>

using namespace TMath;
using namespace std;
namespace rmat {

ClebschGordan::ClebschGordan(double J1, double J2)
{
  /*
  * The constructor needs the angular momenta of the two sub-systems which
  * are coupled. It then calculates the Clebsch-Gordan coefficients and fills 
  * an array that can later be accessed by the user.
  */
  Couple(J1,J2);
}

ClebschGordan::~ClebschGordan(){}

void ClebschGordan::Couple(double J1, double J2)
{
  j1 = J1;
  j2 = J2;
  eps = 1e-4;
  assert(j1 >= 0 && j2 >= 0);
  assert(fmod(2*j1 + eps,1) < 2*eps && fmod(2*j2 + eps,1) < 2*eps);  //Check that the j1,j2 are either integer or half-integer.
  coefficients.clear();
  int Nm1 = (int)lrint(2*j1 + 1);                    //Number of magnetic substates in system 1.
  int Nm2 = (int)lrint(2*j2 + 1);                    //Number of magnetic substates in system 2.
  int NJ = (int)lrint(j1 + j2 - Abs(j1-j2) + 1);     //Number of possible total spins of the combined system.
  coefficients.resize(Nm1);
  for(int i=0; i<Nm1; i++){
    coefficients[i].resize(Nm2);
    for(int j=0; j<Nm2; j++){
      coefficients[i][j].resize(NJ);
      for(int k=0; k<NJ; k++){
        double m1 = -j1 + i;
        double m2 = -j2 + j;
        double J = Abs(j1-j2) + k;
        //coefficients[i][j][k] = CalculateCG(m1,m2,J);
        coefficients[i][j][k] = CalculateCG(j1,j2,J,m1,m2,m1+m2);
        //cout << "(" << m1 << "," << m2 << "|" << J << "," << m1+m2 << ") = " << coefficients[i][j][k] << endl;
      }
    }
  }
}

/*
double ClebschGordan::CGformula(double J1, double J2, double m1, double m2, double J)
{
  double M = m1 + m2;
  assert(lrint(2*J1) >= lrint(2*J2) && lrint(M) >= 0);
  double c = 1.0;
  c *= lrint(2*J + 1);
  c *= Factorial(lrint(J + J1 - J2));
  c *= Factorial(lrint(J - J1 + J2));
  c *= Factorial(lrint(J1 + J2 - J));
  c /= Factorial(lrint(J1 + J2 + J + 1));
  c *= Factorial(lrint(J - M));
  c *= Factorial(lrint(J + M));
  c *= Factorial(lrint(J1 - m1));
  c *= Factorial(lrint(J1 + m1));
  c *= Factorial(lrint(J2 - m2));
  c *= Factorial(lrint(J2 + m2));
  c = Sqrt(c);
  double sum = 0;
  for(int k=0; k<1000; k++){
    if(lrint(J1 + J2 - J - k) < 0 || lrint(J1 -m1 - k) < 0 || lrint(J2 + m2 -k) < 0 ||
       lrint(J - J2 + m1 + k) < 0 || lrint(J - J1 - m2 + k) < 0) continue;
    double term = 1.;
    term *= Power(-1,k);
    term /= Factorial(k);
    term /= Factorial(lrint(J1 + J2 - J - k));
    term /= Factorial(lrint(J1 -m1 - k));
    term /= Factorial(lrint(J2 + m2 -k));
    term /= Factorial(lrint(J - J2 + m1 + k));
    term /= Factorial(lrint(J - J1 - m2 + k));
    sum += term;
  }
  
  c *= sum;
  return c;
}

double ClebschGordan::CalculateCG(double m1, double m2, double J)
{
  if(lrint(m1 + m2) < 0){
    if(lrint(2*j1) < lrint(2*j2)){
      return CGformula(j2,j1,-m2,-m1,J);
    }
    else{
      return CGformula(j1,j2,-m1,-m2,J) * Power(-1,(int)lrint(J - j1 - j2));
    }
  }
  else{
    if(lrint(2*j1) < lrint(2*j2)){
      return CGformula(j2,j1,m2,m1,J) * Power(-1,(int)lrint(J - j1 - j2));
    }
    else{
      return CGformula(j1,j2,m1,m2,J);
    }
  }
}
*/

double ClebschGordan::CalculateCG(double j1, double j2, double j3, double m1, double m2, 
double m3)
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

double ClebschGordan::Coefficient(double m1, double m2, double J)
{
  /*
  * Call this function to extract the calculated Clebsch-Gordan coefficients. If you
  * call it with unphysical arguments it returns zero. The arguments are the magnetic
  * quantum numbers of the two sub-systems and the total angular momentum of the coupled system.
  */

  //First check that we are within proper range.
  if(!((J + eps >= Abs(j1-j2)) && (J - eps <= j1 + j2))) return 0.;
  if(!((m1 + eps >= -j1) && (m1 - eps <= j1))) return 0.;
  if(!((m2 + eps >= -j2) && (m2 - eps <= j2))) return 0.;
  if(!(fmod(2*m1 + eps,1) < 2*eps && fmod(2*m2 + eps,1) < 2*eps && fmod(2*J + eps,1) < 2*eps)) return 0.0;

  int i = lrint(m1 + j1);
  int j = lrint(m2 + j2);
  int k = lrint(J - Abs(j1-j2));

  return coefficients[i][j][k];
}

double ClebschGordan::J1()
{
  return j1;
}

double ClebschGordan::J2()
{
  return j2;
}
}
