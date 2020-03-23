#ifndef SIMULATION_ERROR_H
#define SIMULATION_ERROR_H
#include <memory>
#include <vector>
#include <Math/IFunction.h>
#include <TH1.h>
#include "Simulator.h"

namespace rmat {
namespace threebody {

class SimulationError : public ROOT::Math::IMultiGenFunction {
  /**
  * Given a three-body simulator this class will perform a simulation and 
  * calculate the error with respect to the observed data. An instance of the
  * SimulationError-class can be given to the ROOT::Math::Functor-class and
  * minimised using the Minuit2 minimisation tool. By default it will use the
  * Poisson-likelihood chi-square to quantify the fit error.
  */
  private:
    //Some members have to be mutable, since they are modified in a 'const' fcn.
    mutable Simulator _model;
    std::shared_ptr<TH1> _data;
    mutable int nBins;
    bool normalise;

    virtual double DoEval(const double *x) const override;
    
  public:
    SimulationError(Simulator model, std::shared_ptr<TH1> data);
    ~SimulationError();
    
    //XXX
    virtual ROOT::Math::IBaseFunctionMultiDim *Clone() const override;
    //XXX
    virtual unsigned int NDim() const override;

    virtual double Evaluate() const;
    
    void SetParameters(const double *par) const;
    void PrintParameters();
    
    /**
    * Define whether the simulated spectrum should automatically be normalised 
    * to the number of counts in the observed spectrum. Since the Poisson log-
    * likelihood chi-square is the default statistic and conserves the integral,
    * this option is turned on by default. Can be turned off if another error
    * function is used or the model provides an alternative way of
    * .normalisation
    */
    void AutoNormalise(bool arg = true);
    
    void SetModel(Simulator model);
    void SetData(std::shared_ptr<TH1> data);

    Simulator & GetModel();
    std::shared_ptr<TH1> GetData();
    int Ndf();
    int GetNBins() const;    
};
} //namespace threebody
} //namespace rmat
#endif //SIMULATION_ERROR_H
