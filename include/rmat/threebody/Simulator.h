#ifndef SIMULATOR_H
#define SIMULATOR_H
#include <vector>
#include <memory>
#include <TH1.h>
#include "rmat/threebody/DecayWeight.h"
#include "rmat/threebody/SimEvent.h"

namespace rmat {
namespace threebody {

class Simulator {
  private:
    std::shared_ptr<DecayWeight> weight;
    std::vector<SimEvent> eventPool; //Pool of 'raw' events.
    std::vector<bool> accepted;//Vector indicating which raw events are accepted
    
  public:
    Simulator(std::shared_ptr<DecayWeight> w, std::vector<SimEvent> pool);
    ~Simulator();
    
    /**
    * Run a simulation. Returns the number of accepted events in the simulation.
    */
    long long int Simulate();
    
    /**
    * Obtain the result after the simulation has been performed. Specify the 
    * desired binning of the spectrum.
    */
    std::shared_ptr<TH1> GetSpectrum(int nx, double xMin, double xMax,
                                int ny = 0, double yMin = 0., double yMax = 0.,
                                int nz = 0, double zMin = 0., double zMax = 0.);
                                      
    void SetPool(std::vector<SimEvent> pool);
    void SetDecayWeight(std::shared_ptr<DecayWeight> w);

    std::vector<SimEvent> & GetPool();
    DecayWeight & GetDecayWeight();  

    unsigned int NDim();
    void SetParameters(std::vector<double> par);
};
}//namespace threebody
}//namespace rmat
#endif //SIMULATOR_H
