#ifndef CHEAT_SIMULATOR_H
#define CHEAT_SIMULATOR_H
#include <vector>
#include <memory>
#include <TH1.h>
#include <TH3D.h>
#include "rmat/threebody/DecayWeight.h"
#include "rmat/threebody/Simulator.h"
#include "rmat/threebody/ObservedDalitzPlot.h"

namespace rmat {
namespace threebody {

class CheatSimulator : public Simulator {
  private:
    std::shared_ptr<ObservedDalitzPlot> plot;
    std::shared_ptr<TH1> result; //Filled when running the Simulate() function.

  public:
    CheatSimulator(std::shared_ptr<DecayWeight> w,
                   std::shared_ptr<TH3D> simulated,
                   std::shared_ptr<TH3D> accepted );

    virtual double Simulate();
    
    /**
     * Returns result of calculation. Warning: The dimensions of the returned
     * histogram is determined by the dimensions of the acceptance matrix, and
     * the arguments are ignored.
     */
    virtual std::shared_ptr<TH1> GetSpectrum(int nx, double xMin, double xMax,
                            int ny = 0, double yMin = 0., double yMax = 0.,
                            int nz = 0, double zMin = 0., double zMax = 0.);

};
}//namespace threebody
}//namespace rmat
#endif //CHEAT_SIMULATOR_H
