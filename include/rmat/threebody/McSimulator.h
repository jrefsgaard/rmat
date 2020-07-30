#ifndef MC_SIMULATOR_H
#define MC_SIMULATOR_H
#include <vector>
#include <memory>
#include <TH1.h>
#include "rmat/threebody/DecayWeight.h"
#include "rmat/threebody/Simulator.h"
#include "rmat/threebody/SimEvent.h"

namespace rmat {
namespace threebody {

class McSimulator : public Simulator {
  private:
    std::vector<SimEvent> eventPool; //Pool of 'raw' events.
    std::vector<bool> accepted;//Vector indicating which raw events are accepted
    
  public:
    McSimulator(std::shared_ptr<DecayWeight> w, std::vector<SimEvent> pool);
    ~McSimulator() = default;
    
    virtual double Simulate();
    
    virtual std::shared_ptr<TH1> GetSpectrum(int nx, double xMin, double xMax,
                            int ny = 0, double yMin = 0., double yMax = 0.,
                            int nz = 0, double zMin = 0., double zMax = 0.);
                                      
    void SetPool(std::vector<SimEvent> pool);

    std::vector<SimEvent> & GetPool();
};
}//namespace threebody
}//namespace rmat
#endif //MC_SIMULATOR_H
