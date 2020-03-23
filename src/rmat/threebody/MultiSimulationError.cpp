#include <iostream>
#include <rmat/threebody/MultiSimulationError.h>

using namespace std;
namespace rmat::threebody {

MultiSimulationError::MultiSimulationError() {}

MultiSimulationError::~MultiSimulationError() {}

double MultiSimulationError::DoEval(const double *x) const
{
  SetParameters(x);
  return Evaluate();
}

int MultiSimulationError::AddError(Simulator model, shared_ptr<TH1> data)
{
  if(errors.size() > 0){
    if(model.NDim() != NDim()){
      cout << "MultiSimulationError::AddError(): Mismatch in NDim()!" << endl;
      cout << "  Old NDim() = " << NDim() << ",  requested NDim() = " << model.NDim() << endl;
      exit(EXIT_FAILURE);
    }
  }
  SimulationError err(model,data);
  errors.push_back(err);
  chi2.push_back(0.);
  return errors.size();
}

ROOT::Math::IBaseFunctionMultiDim * MultiSimulationError::Clone() const
{
  MultiSimulationError *clone = new MultiSimulationError(*this);
  return clone;
}
 
unsigned int MultiSimulationError::NDim() const
{
  if(errors.size() > 0){
    return errors.at(0).NDim();
  }
  else return 0;
}

double MultiSimulationError::Evaluate() const
{
  double combined_error = 0.;
  for(int i=0; i<errors.size(); i++){
    double error_i = errors.at(i).Evaluate();
    chi2.at(i) = error_i;
    combined_error += error_i;
    //cout << "MultiSimulationError::Evaluate(): NBins = " << errors.at(i).GetNBins() << endl;
  }
  //cout << "MultiSimulationError::Evaluate(): Error = " << combined_error << endl;
  return combined_error;
}
    
void MultiSimulationError::SetParameters(const double *par) const
{
  for(int i=0; i<errors.size(); i++){
    errors.at(i).SetParameters(par);
  } 
}

int MultiSimulationError::Ndf()
{
  return GetNBins() - NDim();
}

int MultiSimulationError::GetNBins(int index)
{
  if(index == -1){
    int combined_nBins = 0;
    for(int i=0; i<errors.size(); i++){
      combined_nBins += errors.at(i).GetNBins();
    }
    //cout << "MultiSimulationError::GetNBins(): nBins = " << combined_nBins << endl;
    return combined_nBins;
  }
  else{
    return errors.at(index).GetNBins();
  }
}

double MultiSimulationError::GetChi2(int index)
{
  if(index == -1){
    int combined_chi2 = 0;
    for(int i=0; i<chi2.size(); i++){
      combined_chi2 += chi2.at(i);
    }
    //cout << "MultiSimulationError::GetNBins(): nBins = " << combined_nBins << endl;
    return combined_chi2;
  }
  else{
    return chi2.at(index);
  }
}

SimulationError & MultiSimulationError::GetError(int index)
{
  return errors.at(index);
}
}
