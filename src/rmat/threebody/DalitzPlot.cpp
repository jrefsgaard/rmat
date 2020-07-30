#include <rmat/threebody/DalitzPlot.h>
#include <rmat/Nucleus.h>
#include <rmat/RMatrixUtils.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <vector>
#include <iostream>

using namespace std;
using namespace TMath;

namespace rmat::threebody {

DalitzPlot::DalitzPlot(shared_ptr<DecayWeight> w) : weight(w) {}

double DalitzPlot::Value(double x, double y, double Q)
{
  double E23 = 100; //Minimum E23, to cut away 8Be(gs).
  if(y <= -1./Sqrt(3)*x || y >= 1./Sqrt(3)*x || x < 0
    || x*x + y*y >= 0.999
    || y < (Sqrt(3)*x-2.+4.*E23/Q)) return 0.0;

  //if(x*x + y*y >= 0.999) return 0.0;

  //Calculate kinetic energies from Dalitz coordinates.
  double E1 = 0.5 * (x/Sqrt(3) - (y-2)/3) * Q;
  double E2 = (y + 1)/3. * Q;
  double E3 = Q - E1 - E2;

  //cout << "x = " << x << ", y = " << y << ", Q = " << Q << endl;
  //cout << "  E1 = " << E1 << ", E2 = " << E2 << ", E3 = " << E3 << endl;

  double eps1 = 3./2.*E1; //Energy released in first breakup.
  double eps2 = Q - eps1; //Energy released in second breakup.
  
  //cout << "  eps1 = " << eps1 << ", eps2 = " << eps2 << endl;
  double thetap = ACos(Sqrt(3./(eps1*eps2))*(E2 - eps2/2. - eps1/6.)); //Breakup angle in recoil com.
  //cout << "  thetap = " << thetap ;

  //Calculate momenta.
  Nucleus He4(4,2,0,1);
  double alphaMass = He4.M();
  double p1 = Sqrt(2*alphaMass*E1);
  double p2 = Sqrt(2*alphaMass*E2);
  double p3 = Sqrt(2*alphaMass*E3);

  //Create four-vectors with the appropriate length.
  TLorentzVector a1(0,0,p1,E1+alphaMass);
  TLorentzVector a2(0,0,p2,E2+alphaMass);
  TLorentzVector a3(0,0,p3,E3+alphaMass);

  double theta21 = ATan(Sin(thetap)/(Sqrt(eps1/(3*eps2))+Cos(thetap)));
  double theta22 = ATan(Sin(thetap)/(Sqrt(eps1/(3*eps2))-Cos(thetap)));
  //cout << ",  theta21 = " << theta21 << ",  theta22 = " << theta22 << endl;
  a2.RotateY(Pi()-theta21);
  a3.RotateY(Pi()+theta22);
  /*
  TLorentzVector cm = a1 + a2 + a3;
  cout << "  b1 = " << a1.Beta();
  cout << ",  b2 = " << a2.Beta();
  cout << ",  b3 = " << a3.Beta();
  cout << ",  b(cm)= " << cm.Beta();
  cout << ",  K(cm) = " << cm.E() - cm.M() << endl;
  
  cout << "  theta1 = " << a1.Theta() << ",  theta2 = " << a2.Theta();
  cout << ",  theta3 = " << a3.Theta() << endl;
  */
  
  //Create a SimEvent that we can hand over to the weight calculator.
  SimEvent event;
  event.decay.at(0) = a1;
  event.decay.at(1) = a2;
  event.decay.at(2) = a3;
  event.Q = Q;
  event.x = x;
  event.y = y;

  //We calculate the weight for two different event orientations.
  vector<TRotation> v0 {rotations.GetRotation(0)};
  vector<TRotation> v1 {rotations.GetRotation(1)};
  //cout << "Calc. 0:" << endl;
  double w0 = weight->Calculate(event,v0);
  //cout << "Calc. 1:" << endl;
  double w1 = weight->Calculate(event,v1);
  //cout << "w0 = " << w0 << ",  w1 = " << w1 << endl;

  //If the results are different, we average over a hundred different orientations.
  if(IsAlmostZero(w0)) return 0;
  else if(Abs(w1-w0) < 0.01*w0) return w0;
  cout << "DalitzPlot::Value(): Weight not rotation invariant." << endl;
  return weight->Calculate(event,rotations.GetRotations());
}

void DalitzPlot::SetDecayWeight(shared_ptr<DecayWeight> w)
{
  weight = w;
}

DecayWeight & DalitzPlot::GetDecayWeight()
{
  return *(weight.get());
}

unsigned int DalitzPlot::NDim()
{
  return weight->NDim();
}

void DalitzPlot::SetParameters(vector<double> par)
{
  weight->SetParameters(par);
}
}
