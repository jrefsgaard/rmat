#ifndef DALITZ_PLOT_H
#define DALITZ_PLOT_H
#include <memory>
#include "ThreeBodyWeight.h"
#include "../Rotations.h"

namespace rmat {
namespace threebody {

class DalitzPlot {
  private:
    std::shared_ptr<DecayWeight> weight;
    Rotations rotations;

  public:
    DalitzPlot(std::shared_ptr<DecayWeight> w);
    ~DalitzPlot() = default;

    double Value(double x, double y, double Q);

    void SetDecayWeight(std::shared_ptr<DecayWeight> w);

    DecayWeight & GetDecayWeight();  

    unsigned int NDim();
    void SetParameters(std::vector<double> par);
};
}
}
#endif
