#include <iostream>
#include <complex>
#include <cmath>
#include <armadillo>
#include <array>
#include <TLorentzVector.h>
#include <logft/phase_space.h>
#include <rmat/threebody/PeakWeight.h>
#include <rmat/RMatrixUtils.h>

using namespace std;
using namespace arma;
using namespace logft;
namespace rmat::threebody {

PeakWeight::PeakWeight(string type, DoubleCompound s, int peak) : peakID(peak), ThreeBodyWeight(type,s) {}

PeakWeight::~PeakWeight() { }

PeakWeight * PeakWeight::Clone()
{
  return new PeakWeight(*this);
}

void PeakWeight::SetPeakID(int peak)
{
  if(peak < 0 || peak >= scheme.GetNLevels(1)){
    cout << "PeakWeight::SetPeakID(): Peak-index out of range, " << peak << "!" << endl;
    exit(EXIT_FAILURE);
  }
  peakID = peak;
}

int PeakWeight::GetPeakID()
{
  return peakID;
}

double PeakWeight::Calculate(SimEvent &event)
{
  //Return 0 if this observed Q-value has been excluded from the calculation.
  double Qobs = event.Q;
  for(auto range : excludedQs){
    if(range[0] < Qobs && Qobs < range[1]){
      return 0.;
    }
  }

  array<TLorentzVector,3> &decay = event.decay;

  double Q = 0; //True Q-value of the alpha breakup.
  for(auto & alpha : decay){ Q += (alpha.Energy() - alpha.M());}

  //Calculate excitation energy of primary system.
  double Ex1 = Q + thresholds.at(0);
  
  double El2 = scheme.GetEnergy(1,peakID); //Excitation energy of level b.
  double Eth2 = thresholds.at(1);
  double Ec = Q - (El2 - Eth2); //Channel energy in primary channel.
  
  //cout << "Ex1 = " << Ex1 << ",  Eth2 = " << Eth2 << endl;
  //cout << "Q = " << Q << ",  Ec = " << Ec << endl;

  //We determine the maximum Ja we need to include in the sum.
  int JaMax = ReturnMaxSpin(scheme.GetJs(0));
  //cout << "JaMax = " << JaMax << endl;
   
  //Mat<complex<double>> amplitudes(JaMax+1,2*JaMax+1);
  double weight = 0.;

  //cout << "line 171" << endl;
  
  //Then we start the proper calculation.
  for(int Ja = 0; Ja <= JaMax; Ja++){
    //cout << "line 175" << endl;
    Mat<complex<double>> A1 = scheme.LevelMatrix(0,Ex1,Ja);  //Primary level matrix.
    //cout << "line 177" << endl;
    if(A1.size() == 0) continue; //Skip if A1 is empty (for instance if no levels of spin Ja exist).
    //cout << "line 179" << endl;
    //cout << "A1 = " << A1 << endl;

    //The indices of the levels with the proper spin.
    vector<int> pIndices = scheme.GetLevelIndices(0,Ja);
        
    //Sum over l1
    int Nl1 = scheme.GetPrimaryChannels().size();
    //cout << "Nl1 = " << Nl1 << endl;
    for(int il=0; il<Nl1; il++){
      Channel &ci = scheme.GetPrimaryChannels().at(il);
      double P1 = ci.Penetrability(Ec);
      //cout << "P1 = " << P1 << endl;
      int c = scheme.GetChannelIndex(il,peakID);
      
      complex<double> amplitude(0.,0.);
      //cout << "il = " << il << ",  c = " << c << endl;
      //We use the index vectors to only loop over levels with the proper J.
      for(int m=0; m<pIndices.size(); m++){    //Level index with respect to other levels of the same spin.
        int lambda = pIndices.at(m);           //Global level index
        double Blambda = Bbeta.at(lambda);
        if(IsAlmostZero(Blambda)) continue;
        //cout << "m = " << m << ",  lambda = " << lambda << endl;
        for(int n=0; n<pIndices.size(); n++){
          int mu = pIndices.at(n);
          //cout << "lambda = " << lambda << ",  mu = " << mu << endl;
          double gammaMu = scheme.GetWidthAmplitude(0,mu,c);
          //cout << "gammaLambda = " << gammaLambda << endl;
          if(IsAlmostZero(gammaMu)) continue;
          complex<double> f = A1(m,n) * Blambda * gammaMu;
          amplitude += f;
        }
      }
      weight += P1 * norm(amplitude);
      //cout << "Weight = " << P1 * norm(amplitude) << endl;
    }
  }
  
  //Calculate the beta-decay phase-space factor.
  double fBeta = 0.;
  int A = scheme.GetA(0);
  int Zf = scheme.GetZ(0);
  int Zi = Zf - dZ;
  try { fBeta = logft::calculatePhaseSpace(Zi,Zf,A,Ex1);} 
  catch (...) { fBeta = 0.;}             
  //fBeta = 1.;
  return fBeta * weight;
}
} //namespace rmat::threebody;
