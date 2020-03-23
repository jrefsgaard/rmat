#include <iostream>
#include <TMath.h>
#include <rmat/threebody/CoulombCorrector.h>

using namespace std;

namespace rmat::threebody {

CoulombCorrector::CoulombCorrector() {}

CoulombCorrector::CoulombCorrector(DoubleCompound &scheme, double radius)
{
  for(Channel c : scheme.GetPrimaryChannels()){
    //cout << "primary channel, L = " << c.L() << endl;
    c.SetRadius(radius);
    pChannels.push_back(c);
  }
  for(Channel c : scheme.GetSecondaryChannels()){
    //cout << "secondary channel, L = " << c.L() << endl;
    c.SetRadius(radius);
    sChannels.push_back(c);
  }
}

CoulombCorrector::~CoulombCorrector() {}

double CoulombCorrector::GetFSCI(double E, double E23, double E12, double E13, int pIndex, int sIndex)
{
  if(pChannels.size() == 0 || sChannels.size() == 0) return 1.;

  Channel &cp = pChannels.at(pIndex);
  Channel &cs = sChannels.at(sIndex);

  double correction = 1.;
  correction *= cp.Rho(E) / cp.Penetrability(E);
  correction *= cs.Penetrability(E12) / cs.Rho(E12);
  correction *= cs.Penetrability(E13) / cs.Rho(E13);

  /*
  static int n = 0;
  n++;
  if(TMath::Sqrt(correction) > 10.){
    cout << "FSCI = " << correction << ", n = " << n << ",  E1-23 = " << E << ",  Etot = " << E + E23 << endl;
    n = 0;
    cout << "rho1 = " << cp.Rho(E) << ",  P1 = " << cp.Penetrability(E) << ",  rho12 = " << cs.Rho(E12) << ",  P12 = " << cs.Penetrability(E12);
    cout <<  ",  rho13 = " << cs.Rho(E13) << ",  P13 = " << cs.Penetrability(E13) << endl;
  }
  */
  return correction;
}

}
