#ifndef WEIGHT_CALCULATOR_H
#define WEIGHT_CALCULATOR_H
#include <memory>
#include <vector>
#include <array>
#include <TLorentzVector.h>
#include "rmat/threebody/DecayWeight.h"
#include "rmat/threebody/SimEvent.h"

namespace rmat::threebody {

class WeightCalculator {
  private:
    std::shared_ptr<DecayWeight> weight;
    double maxWeight;
    
  public:
    WeightCalculator(std::shared_ptr<DecayWeight> w);
    ~WeightCalculator();
    
    /**
    * This function is potentially very heavy computationally, and it is there-
    * fore employing OpenMP for parrallelisation.
    */
    std::vector<double> CalculateWeights(std::vector<SimEvent> & raw);
    
    //Get maximum weight from last run of CalculateWeights().
    double MaxWeight();
    
    void SetWeight(std::shared_ptr<DecayWeight> w);
    DecayWeight & GetWeight();

};
}//namespace rmat::threebody
#endif //WEIGHT_CALCULATOR_H
