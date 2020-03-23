#ifndef SAMPLER_H
#define SAMPLER_H
#include <vector>
#include <TRandom3.h>
#include <rmat/threebody/SimEvent.h>

namespace rmat {
namespace threebody {

class Sampler {
  private:
    TRandom3 rGen;
    
  public:
    Sampler();
    ~Sampler();
    
    /**
    * Returns a vector with accepted events.
    */
    std::vector<SimEvent> Sample(std::vector<SimEvent> & raw, std::vector<double> & weights, double maxWeight);
    
    /**
    * Returns a vector with boolean values specifying whether a given event was
    * accepted or not. Returned vector has the same length as the weight vector.
    */
    std::vector<bool> Sample(std::vector<double> & weights, double maxWeight);
    
    long long int Sample(std::vector<double> & weights, double maxWeight, std::vector<bool> & accepted);

    /**
    * A simple function which will accept or reject an event based on the given
    * weight. It assumes uniform weight distribution between 0 and maxWeight.
    */    
    bool Accept(double weight, double maxWeight);
};
} //namespace threebody;
} //namespace rmat;
#endif //SAMPLER_H
