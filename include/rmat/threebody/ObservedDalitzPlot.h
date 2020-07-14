#ifndef OBSERVED_DALITZ_PLOT_H
#define OBSERVED_DALITZ_PLOT_H
#include "DalitzPlot.h"
#include <memory>
#include <vector>
#include <TH3D.h>

namespace rmat {
namespace threebody {

class ObservedDalitzPlot {
  private:
    std::shared_ptr<DalitzPlot> model;
    std::shared_ptr<TH3D> acceptance;

  public:
    ObservedDalitzPlot(std::shared_ptr<DalitzPlot> model_,
                       std::shared_ptr<TH3D> simulated,
                       std::shared_ptr<TH3D> accepted);
    ~ObservedDalitzPlot() = default;

    void SetModel(std::shared_ptr<DalitzPlot> model_);
    DalitzPlot & GetModel();

    std::shared_ptr<TH3D> GetSpectrum();

    void SetParameters(std::vector<double> par);
    void PrintParameters();
    int NDim();
};
}
}
#endif
