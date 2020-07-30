#ifndef MULTI_SIMULATION_ERROR_H
#define MULTI_SIMULATION_ERROR_H
#include <vector>
#include <memory>
#include <Math/IFunction.h>
#include <TH1.h>
#include "Simulator.h"
#include "SimulationError.h"

namespace rmat {
namespace threebody {

class MultiSimulationError : public ROOT::Math::IMultiGenFunction {
  /**
  * Class to hold several SimulationErrors. Behaves otherwise exactly
  * as a normal SimulationError.
  */
  private:
    std::vector<SimulationError> errors;
    mutable std::vector<double> chi2;
    virtual double DoEval(const double *x) const override;

  public:
    MultiSimulationError();
    ~MultiSimulationError();

    int AddError(std::shared_ptr<Simulator> model, std::shared_ptr<TH1> data);
    SimulationError & GetError(int index);

    //XXX
    virtual ROOT::Math::IBaseFunctionMultiDim *Clone() const override;
    //XXX
    virtual unsigned int NDim() const override;

    virtual double Evaluate() const;
    
    void SetParameters(const double *par) const;
    int Ndf();
    int GetNBins(int index = -1);
    double GetChi2(int index = -1);  
};
} //namespace threebody
} //namespace rmat
#endif //MULTI_SIMULATION_ERROR_H
