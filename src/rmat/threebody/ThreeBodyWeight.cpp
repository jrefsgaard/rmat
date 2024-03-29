#include <iostream>
#include <complex>
#include <cmath>
#include <armadillo>
#include <array>
#include <TLorentzVector.h>
//#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <logft/phase_space.h>
#include <rmat/threebody/ThreeBodyWeight.h>
#include <rmat/RMatrixUtils.h>

using namespace std;
//using namespace boost::math;
using namespace arma;
using namespace logft;
namespace rmat::threebody {

ThreeBodyWeight::ThreeBodyWeight(std::string type, DoubleCompound s) : scheme(s)
{
  SetType(type);
  Bbeta.resize(scheme.GetNLevels(0),0.);
  thresholds.resize(2);
  doCorrection = false;
  //Construct the angular momentum couplings.
  vector<Channel> & channels = scheme.GetPrimaryChannels();
  vector<Level> & levels = scheme.GetSystem(1).GetLevels();
  cgCoefficients.resize(channels.size());
  int Lmax = -1;
  for(int c=0; c<channels.size(); c++){
    cgCoefficients.at(c).resize(levels.size());
    for(int mu=0; mu<levels.size(); mu++){
      int L1 = channels.at(c).L();
      double Jb = levels.at(mu).J();
      ClebschGordan cg(L1,Jb);
      cgCoefficients.at(c).at(mu) = cg;
      if(L1 > Lmax) Lmax = L1;
    }
  }
}

ThreeBodyWeight::~ThreeBodyWeight() { }

//Std. ordering: First parameters for primary compound system, then for secondary system.
//{{E1,theta11,theta12,..,B1,E2,theta21,...},{E1,Gamma11,Gamma12,....,E2,Gamma21,...}}
void ThreeBodyWeight::SetParameters(vector<double> par)
{
  if(par.size() != NDim()){
    cout << "  ThreeBodyWeight::SetParameters(): Wrong number of parameters, par.size() = "<< par.size() << "!" << endl;
    exit(EXIT_FAILURE);
  }
  //First, we set the parameters for the primary compound system.
  int Nla = scheme.GetNLevels(0);
  int Nc1 = scheme.GetNChannels(0);
  for(int l=0; l<Nla; l++){
    double El = par.at(l*(Nc1 + 2));
    double Bl = par.at((l+1)*(Nc1+2)-1);
    /*
    One could use partial widths here, but we prefer dimensionless widths.
    vector<double> Gammal;
    for(int c=0; c<Nc1; c++){
      double Gammalc = par.at(l*(Nc1+2) + 1 + c);
      Gammal.push_back(Gammalc);
    }
    //cout << El << "  " << Gammal.at(0) << "  " << Gammal.at(1) << "  " << Bl << endl;
    */
    vector<double> thetal; //vector of theta^2
    for(int c=0; c<Nc1; c++){
      double thetalc = par.at(l*(Nc1+2) + 1 + c);
      thetal.push_back(thetalc);
    }
    scheme.SetEnergy(0,l,El);
    scheme.SetDimensionlessWidths(0,l,thetal);
    //scheme.SetPartialWidths(0,l,Gammal);
    SetB(l,Bl);
  }
  int N1 = Nla * (Nc1 + 2); //The number of parameters describing the primary system.
  //Then the secondary compound system.
  int Nlb = scheme.GetNLevels(1);
  int Nc2 = scheme.GetNChannels(1);
  for(int l=0; l<Nlb; l++){
    double El = par.at(N1 + l*(Nc2 + 1));
    vector<double> Gammal;
    for(int c=0; c<Nc2; c++){
      double Gammalc = par.at(N1 + l*(Nc2 + 1) + 1 + c);
      Gammal.push_back(Gammalc);
    }
    //cout << El << "  " << Gammal.at(0) <<  endl;
    scheme.SetEnergy(1,l,El);
    scheme.SetPartialWidths(1,l,Gammal);
  }
}

void ThreeBodyWeight::PrintParameters()
{
  cout << endl;
  cout << "ThreeBodyWeight: Printing on-resonance R-matrix parameters..." << endl;
  vector<string> systems {"Primary system:", "Secondary system:"};
  for(int id=0; id<=1; id++){
    cout << systems.at(id) << endl;
    int Nl = scheme.GetNLevels(id);
    int Nc = scheme.GetNChannels(id);
    for(int l=0; l<Nl; l++){
      cout << "  Level " << l << ":" << endl;
      cout << "    E = " << scheme.GetEnergy(id,l) << endl;
      cout << "    Gamma = " << scheme.GetWidth(id,l) << endl;
      for(int c=0; c<Nc; c++){
        cout << "      g(" << l << c << ") = " << scheme.GetWidthAmplitude(id,l,c) << endl;
        cout << "      Gamma(" << l << c << ") = " << scheme.GetPartialWidth(id,l,c) << endl;
        cout << "      theta(" << l << c << ") = " << scheme.GetDimensionlessWidth(id,l,c) << endl;
        //cout << "      Wigner limit = " << scheme.GetSystem(id).GetChannel(c).WignerLimit() << endl;
      }
      if(id == 0) cout << "    B(" << l << ") = " << Bbeta.at(l) << endl;
    }
  }
  cout << endl ;
}

int ThreeBodyWeight::NDim()
{
  int Nla = scheme.GetNLevels(0);
  int Nlb = scheme.GetNLevels(1);
  int Nc1 = scheme.GetNChannels(0);
  int Nc2 = scheme.GetNChannels(1);
  return Nla * (Nc1 + 2) + Nlb * (Nc2 +1);
}
    
void ThreeBodyWeight::SetType(std::string type)
{
  if(type == "-") dZ = 1;
  else if(type == "+") dZ = -1;
  else{
    cout << "  ThreeBodyWeight::SetType(): Unknown decay type: "<< type << "!" << endl;
    exit(EXIT_FAILURE);
  }
}

void ThreeBodyWeight::EnsureSize(std::vector<double> & B)
{
  if(B.size() != scheme.GetNLevels(0)){
    cout << "  ThreeBodyWeight::EnsureSize(): Wrong number of feeding factors. Resizing..." << endl;
    B.resize(scheme.GetNLevels(0),0.);
  }
}

void ThreeBodyWeight::SetBs(std::vector<double> B)
{
  EnsureSize(B);
  Bbeta = B;
}

void ThreeBodyWeight::SetB(int lambda, double B)
{
  if(lambda < 0 || lambda >= Bbeta.size()){
    cout << "  ThreeBodyWeight::SetB(): Index out of range: "<< lambda << "!" << endl;
    exit(EXIT_FAILURE);
  }
  Bbeta.at(lambda) = B;
}

DoubleCompound & ThreeBodyWeight::GetScheme()
{
  return scheme;
}

void ThreeBodyWeight::SetThresholds(double E1, double E2)
{
  thresholds.at(0) = E1;
  thresholds.at(1) = E2;
}

vector<double> & ThreeBodyWeight::GetThresholds()
{
  return thresholds;
}

int ThreeBodyWeight::ReturnMaxSpin(std::vector<double> spins)
{
  int JMax = -1.;
  for(double Ji : spins){
    int Jx2 = lrint(2*Ji);
    if(Jx2 % 2 != 0){
      cout << "  ThreeBodyWeight::ReturnMaxSpin(): Calculation only handles integer spins, J: "<< Ji << "!" << endl;
      exit(EXIT_FAILURE);
    }
    int Jint = lrint(Ji);
    JMax = Jint > JMax ? Jint : JMax;
  }
  return JMax;
}

ThreeBodyWeight * ThreeBodyWeight::Clone()
{
  return new ThreeBodyWeight(*this);
}

void ThreeBodyWeight::DoCoulombCorrection(double radius)
{
  CoulombCorrector newCorrector(scheme,radius);
  coulombCorrector = newCorrector;
}
/*
double ThreeBodyWeight::Calculate(SimEvent &event)
{
  vector<TRotation> v {TRotation()}; //Identity rotation.
  return Calculate(event,v);
}
*/
double ThreeBodyWeight::Calculate(SimEvent &event)
{
  //Return 0 if this observed Q-value has been excluded from the calculation.
  double Qobs = event.Q;
  for(auto range : excludedQs){
    if(range[0] < Qobs && Qobs < range[1]){
        //cout << "Q-value = " << Q << endl;
      return 0.;
    }
  }

  array<TLorentzVector,3> &decay = event.decay;
  
  //True Q-value of the alpha breakup.
  double Q = 0;
  for(auto & alpha : decay){ Q += (alpha.Energy() - alpha.M());}
  
  //Calculate excitation energy of primary system.
  double Ex1 = Q + thresholds.at(0);

  //Also check that all spins are integer.
  int JaMax = ReturnMaxSpin(scheme.GetJs(0));
 
  //Same for Jb
  int JbMax = ReturnMaxSpin(scheme.GetJs(1));
  
  Mat<complex<double>> amplitudes(JaMax+1,2*JaMax+1);
  complex<double> zero(0.,0.);
  complex<double> i(0.,1.);
  amplitudes.fill(zero);
  
  //Then we start the proper calculation.
  for(int Ja = 0; Ja <= JaMax; Ja++){
    //Primary level matrix.
    Mat<complex<double>> A1 = scheme.LevelMatrix(0,Ex1,Ja); 
    //Skip if A1 is empty (for instance if no levels of spin Ja exist).
    if(A1.size() == 0) continue;

    //The indices of the levels with the proper spin.
    vector<int> pIndices = scheme.GetLevelIndices(0,Ja);
      
    //Symmetrisation in the alpha particle labels.
    for(int j=0; j<3; j++){
      TLorentzVector alpha1 = decay.at(j);
      TLorentzVector alpha2 = decay.at((j+1)%3);
      TLorentzVector alpha3 = decay.at((j+2)%3);  
    
      //We find angles and energies.
      double E1 = (alpha1.Energy() - alpha1.M());
      double theta1angle = alpha1.Theta();
      double phi1angle = alpha1.Phi();  //Mind confusion with HS phase shift!

      TLorentzVector rcm = alpha2 + alpha3;  //Recoil center of mass system.
      alpha2.Boost(-rcm.BoostVector());
      double theta2 = alpha2.Theta();
      double phi2 = alpha2.Phi();

      rcm.Boost(-rcm.BoostVector());
      double E23 = (rcm.Energy() - alpha2.M() - alpha3.M());//Relative energy of alpha2 and alpha3.
      double Ex2 = E23 + thresholds.at(1);
      double E = Q - E23;  //Relative energy of alpha1 and recoil.
     
      //For the Coulomb correction we need the relative energies of alpha 1-2 and alpha 1-3.
      TLorentzVector cm12 = alpha1 + alpha2;
      cm12.Boost(-cm12.BoostVector());
      double E12 = cm12.Energy() - alpha1.M() - alpha2.M();
      TLorentzVector cm13 = alpha1 + alpha3;
      cm13.Boost(-cm13.BoostVector());
      double E13 = cm13.Energy() - alpha1.M() - alpha3.M(); 
 
      for(int Jb=0; Jb <= JbMax; Jb++){
        Mat<complex<double>> A2 = scheme.LevelMatrix(1,Ex2,Jb);
        if(A2.size() == 0) continue; //Skip if matrix is empty.
        //A2(0,0) = complex<double>(1,0);//XXX

        //Secondary levels with proper spin.
        vector<int> sIndices = scheme.GetLevelIndices(1,Jb);
        
        //Sum over l1
        int Nl1 = scheme.GetPrimaryChannels().size();
        for(int il=0; il<Nl1; il++){
          Channel &ci = scheme.GetPrimaryChannels().at(il);
          int L1 = ci.L();
          double P1 = ci.Penetrability(E);
          double rho1 = ci.Rho(E);
          double phi1 = ci.HardSphereShift(E);
          double omega1 = ci.CoulombShift(E);
          complex<double> Omega1 = exp(i * (omega1 - phi1));
                
          //Sum over L2
          int Nl2 = scheme.GetSecondaryChannels().size();
          for(int ilp=0; ilp<Nl2; ilp++){
            double FSCI = coulombCorrector.GetFSCI(E,E23,E12,E13,il,ilp);
            Channel &cip = scheme.GetSecondaryChannels().at(ilp);
            int L2 = cip.L();
            double P23 = cip.Penetrability(E23);
            double rho23 = cip.Rho(E23);
            double phi23 = cip.HardSphereShift(E23);
            double omega23 = cip.CoulombShift(E23);
            complex<double> Omega23 = exp(i * (omega23 - phi23));
            
            //We use the index vectors to only loop over levels with the proper J.
            for(int m=0; m<pIndices.size(); m++){    //Level index with respect to other levels of the same spin.
              int lambda = pIndices.at(m);           //Global level index
              //cout << "m = " << m << ",  lambda = " << lambda << endl;
              for(int n=0; n<pIndices.size(); n++){
                int mu = pIndices.at(n);
                //cout << "n = " << n << ",  mu = " << mu << endl;
                double Bmu = Bbeta.at(mu);
                if(IsAlmostZero(Bmu)) continue;
                for(int mp=0; mp<sIndices.size(); mp++){
                  int lambdap = sIndices.at(mp);
                  //cout << "mp = " << mp << ",  lambdap = " << lambdap << endl;
                  double gammaLambdap = scheme.GetWidthAmplitude(1,lambdap,ilp);
                  //cout << "gammaLambdap = " << gammaLambdap << endl;
                  if(IsAlmostZero(gammaLambdap)) continue;
                  for(int np=0; np<sIndices.size(); np++){  
                    int mup = sIndices.at(np);
                    //cout << "np = " << np << ",  mup = " << mup << endl;
                    int c = scheme.GetChannelIndex(il,mup);
                    ClebschGordan &cg = cgCoefficients.at(il).at(mup);

                    double gammaLambda = scheme.GetWidthAmplitude(0,lambda,c);
                    //cout << "gammaLambda = " << gammaLambda << endl;
                    if(IsAlmostZero(gammaLambda)) continue;
                    for(int ma=-Ja; ma<=Ja; ma++){
                      //cout << ",  ma = " << ma ;
                      for(int mb = -Jb; mb <= Jb; mb++){
                        if(abs(ma-mb) > L1) continue;
                        
                        //Angular part.
                        complex<double> Cmmj(cg.Coefficient(ma-mb,mb,Ja),0);
                        complex<double> Ylm1(sphericalHarmonic.Value(L1,ma-mb,theta1angle,phi1angle));
                        complex<double> Ylm2(sphericalHarmonic.Value(L2,mb,theta2,phi2));
                        complex<double> f = Cmmj * Ylm1 * Ylm2 / sqrt(2*Ja+1);
                        //f *= pow(i,L1+L2); //Good for T-inversal invariance..?
                        
                        //Primary breakup
                        f *= Omega1 * sqrt(P1);
                        
                        //Intermediate resonance part
                        f *= Omega23 * sqrt(P23) * gammaLambdap * A2(mp,np);
                        f /= sqrt(M_PI);
                        
                        //Coordinate transformation part.
                        f /= pow(E*E23,0.25);
                        
                        //Normalisation factor for symmetrisation.
                        f /= sqrt(3);
                        
                        //Properties of 12C resonance.
                        f *= Bmu * A1(m,n) * gammaLambda;
                        
                        //Ad-hoc final-state Coulomb interaction correction
                        f *= sqrt(FSCI);
                        
                        amplitudes(Ja,ma+JaMax) += f;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  double weight = sum(sum(square(abs(amplitudes))));
  
  //Calculate the beta-decay phase-space factor.
  double fBeta = 0.;
  int A = scheme.GetA(0);
  int Zf = scheme.GetZ(0);
  int Zi = Zf - dZ;
  try { fBeta = logft::calculatePhaseSpace(Zi,Zf,A,Ex1);} 
  catch (...) { fBeta = 0.;}             
  //fBeta = 1.;

  double PS = 2 * pow(M_PI*Q,2);
  return fBeta * PS * weight;
}
} //namespace rmat::threebody;
