#ifndef THREE_BODY_WEIGHT_H
#define THREE_BODY_WEIGHT_H
#include <vector>
//#include <array>
#include <string>
//#include <TLorentzVector.h>
#include "DecayWeight.h"
#include "DoubleCompound.h"
#include "CoulombCorrector.h"
#include "SimEvent.h"
#include "../ClebschGordan.h"
#include "../SphericalHarmonic.h"
#include <TRotation.h>

namespace rmat {
namespace threebody {

class ThreeBodyWeight : public DecayWeight {
  protected:
    int dZ;
    std::vector<double> thresholds; //3body threshold
    std::vector<double> Bbeta;
    DoubleCompound scheme;
    bool doCorrection;
    CoulombCorrector coulombCorrector;    
    std::vector<std::vector<ClebschGordan>> cgCoefficients;
    SphericalHarmonic sphericalHarmonic;

    void EnsureSize(std::vector<double> & B);  //Ensures that the vector has proper length.
    int ReturnMaxSpin(std::vector<double> spins);
    
  public:
    ThreeBodyWeight(std::string type, DoubleCompound s);
    ~ThreeBodyWeight();
    
    //virtual double Calculate(std::array<TLorentzVector,3> &);
    //virtual double Calculate(SimEvent &event);
    virtual double Calculate(SimEvent &event, std::vector<TRotation> &rotations);
    
    /**
    * Std. ordering: First parameters for primary compound system, then for secondary system.
    * {{E1,theta11,theta12,..,B1,E2,theta21,...},{E1,Gamma11,Gamma12,....,E2,Gamma21,...}}
    * Note: The theta's are actually theta^2 and the ordering is:
    * {(l,lambda),(l,lambda'),...,(l',lambda),(l',lambda'),...}.
    */
    virtual void SetParameters(std::vector<double>);
    virtual void PrintParameters();
    int NDim();
    
    void SetType(std::string type);
    void SetBs(std::vector<double> B);
    void SetB(int lambda, double B);
    void SetThresholds(double E1, double E2);

    void DoCoulombCorrection(double radius);
    
    DoubleCompound & GetScheme();
    std::vector<double> & GetThresholds();
    
    virtual ThreeBodyWeight * Clone();
};
} //namespace rmat;
}//naespace threebody;
#endif //THREE_BODY_WEIGHT_H
