#ifndef MULTI_ERROR_H
#define MULTI_ERROR_H
#include <Math/IFunction.h>
#include <vector>
#include <memory>

namespace rmat {

/**
 * This class can hold several functions that need to be minimised
 * simultaneously.
 */
class MultiError : public ROOT::Math::IMultiGenFunction {
  private:
    std::vector<std::shared_ptr<ROOT::Math::IMultiGenFunction>> errors;

    //Mandatory for IMultiGenFunction
    virtual double DoEval(const double *x) const override;

  public:
    MultiError() = default;
    ~MultiError() = default;

    void AddError(std::shared_ptr<ROOT::Math::IMultiGenFunction> error);

    ROOT::Math::IBaseFunctionMultiDim & GetError(int i);

    //Mandatory for IMultiGenFunction
    virtual ROOT::Math::IBaseFunctionMultiDim *Clone() const override;

    //Mandatory for IMultiGenFunction.
    //Returns dimension of the error with the largest dimension.
    virtual unsigned int NDim() const override;
};
}
#endif
