#ifndef PHASE_SPACE_GENERATOR_H
#define PHASE_SPACE GENERATOR_H
#include <vector>
#include <TRandom3.h>
#include "rmat/threebody/SimEvent.h"
#include "rmat/threebody/EventGenerator.h"
#include "rmat/threebody/Sampler.h"

namespace rmat {
namespace threebody {

class PhaseSpaceGenerator {
  private:
    TRandom3 rGen;
    EventGenerator generator;
    Sampler sampler;
    
  public:
    PhaseSpaceGenerator();
    ~PhaseSpaceGenerator();
    
    std::vector<SimEvent> Generate(int N, double Q);
    
    std::vector<SimEvent> Generate(int N, double Qmin, double Qmax);

};
} //namespace threebody;
} //namespace rmat;
#endif //PHASE_SPACE_GENERATOR_H
