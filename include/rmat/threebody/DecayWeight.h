#ifndef DECAY_WEIGHT_H
#define DECAY_WEIGHT_H
#include <vector>
#include <array>
#include "SimEvent.h"

namespace rmat {
namespace threebody {

class DecayWeight {
  protected:
    std::vector<std::array<double,2>> excludedQs; //These Q(3alpha) ranges should be excluded.

  public:
    DecayWeight() = default;
    virtual ~DecayWeight() = default;

    //virtual double Calculate(std::array<TLorentzVector,3> &) = 0;
    virtual double Calculate(SimEvent &) = 0;
    virtual void SetParameters(std::vector<double>) = 0;
    virtual void PrintParameters() = 0;
    virtual int NDim() = 0;
    
    virtual DecayWeight * Clone() = 0;

    /**
    * It is the responsibility of the child class to implement the actual exclusion
    * in the calculation of the decay weight.
    */    
    void ExcludeQRange(double Qmin, double Qmax);
    std::vector<std::array<double,2>> & GetExcludedQRanges();
};
} //namespace threebody;
} //namespace rmat;
#endif //DECAY_WEIGHT_H
