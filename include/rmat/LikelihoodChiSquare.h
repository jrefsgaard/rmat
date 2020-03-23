#ifndef LIKELIHOOD_CHI_SQUARE_H
#define LIKELIHOOD_CHI_SQUARE_H
#include <memory>
#include <Math/IFunction.h>
#include "Spectrum.h"
#include "ChiSquare.h"

namespace rmat {

class LikelihoodChiSquare : public ChiSquare {
  public:
    LikelihoodChiSquare(std::shared_ptr<TH1> d, std::shared_ptr<Spectrum> m);
    ~LikelihoodChiSquare();
    
    virtual ROOT::Math::IBaseFunctionMultiDim *Clone() const override;
    virtual double Evaluate() const;
};
} //namespace rmat;
#endif //LIKELIHOOD_CHI_SQUARE_H
