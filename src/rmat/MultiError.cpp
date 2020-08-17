#include <rmat/MultiError.h>

namespace rmat {

void MultiError::AddError(std::shared_ptr<ROOT::Math::IMultiGenFunction> error)
{
  errors.push_back(error);
}

ROOT::Math::IBaseFunctionMultiDim & MultiError::GetError(int i)
{
  return *(errors.at(i).get());
}

double MultiError::DoEval(const double *x) const
{
  double chi = 0;
  for(auto error : errors) chi += (*error.get())(x);
  return chi;
}

ROOT::Math::IBaseFunctionMultiDim *MultiError::Clone() const
{
  MultiError *clone = new MultiError(*this);  
  return clone;
}

unsigned int MultiError::NDim() const
{
  unsigned int dim = 0;
  for(auto error : errors) dim = error->NDim() > dim ? error->NDim() : dim;
  return dim; 
}
}
