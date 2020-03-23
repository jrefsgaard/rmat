#include "rmat/CompoundSystem.h"

namespace rmat {
namespace threebody {

class DensityOfStates {
  private:
    CompoundSystem system;
    
  public:
    DensityOfStates(int A, int Z);
    ~DensityOfStates();
    
    //Observed level energy and partial width.
    void SetLevel(double E, double Gamma);

    void SetChannel(Channel c, double threshold);
    void SetChannel(std::string type, int l, double r0, double threshold);
    
    void UseInterpolation(bool interp = true, double Emax = 10000., double Estep = 20.);
    double Evaluate(double Ec);
};
}//namespace threebody
}//namespace rmat
