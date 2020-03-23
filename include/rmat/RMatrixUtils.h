#ifndef RMATRIX_UTILS_H
#define RMATRIX_UTILS_H
#include <complex>

namespace rmat {
//Main functions.
double P(int, double, double);
double ShiftNeg(int, double, double);
double ShiftPos(int, double, double);

//Helper functions.
double CoulombF(int, double, double);
double CoulombG(int, double, double);
int CoulombFunctions(int, double, double, double&, double&);
int CoulombFunctions(int, double, double, double&, double&, double&, double&);
double WhittakerW(double, double, double);
struct Wparams;
double fWhittakerW(double, void *);
double WhittakerWDiff(double, double, double);
bool CheckFinity(std::complex<double>&, std::complex<double>&, std::complex<double>&, std::complex<double>&);
int Delta(int, int);
bool IsAlmostZero(double num, double smudgeFactor = 10.);
double IsAlmostEqual(double num1, double num2);
double ClebschGordanFCN(double j1, double j2, double j3, double m1, double m2, double m3);
std::complex<double> SphericalHarmonic(int l, int m, double theta, double phi);
} //namespace rmat
#endif //RMATRIX_UTILS
