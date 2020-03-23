#ifndef METROPOLIS_WALKER_H
#define METROPOLIS_WALKER_H
#include <vector>
#include "Math/Functor.h"

namespace rmat {
struct Step {
  std::vector<double> x;
  double value;
};

struct ParameterInfo {
  bool fixed;
  double sigma;
  double min;
  double max;
};

class MetropolisWalker {
  private:
    ROOT::Math::Functor chi2;            //The function to sample.
    double T;                            //Temperature of the chi2 function.
    int Nscale;                          //Steps between scaling of sigmas.
    std::vector<Step> position;          //All accepted steps.
    std::vector<Step> history;           //All evaluations.
    std::vector<ParameterInfo> options;  //Options for sampling the parameters.
    
    double Evaluate(std::vector<double> x);   //Evaluate and store result in history.
    
  public:
    MetropolisWalker(ROOT::Math::Functor f, int scaling = 10);
    ~MetropolisWalker();

    void Walk(std::vector<double> x0, int steps);
    
    void SetFunction(ROOT::Math::Functor f);
    void SetTemperature(double temp);
    void SetScalingPeriod(int nsteps);
    void SetParOptions(int ipar, double sigma, double pmin, double pmax);
    void SetParOptions(std::vector<ParameterInfo> opts);
    void FixParameter(int ipar);
    void ReleasePar(int ipar);
    
    ROOT::Math::Functor & GetFunction();
    double GetTemperature();
    int GetScalingPeriod();
    std::vector<Step> & GetSteps();
    std::vector<Step> & GetHistory();
    ParameterInfo & GetParOptions(int ipar);
    std::vector<ParameterInfo> & GetParOptions();
};
}
#endif
