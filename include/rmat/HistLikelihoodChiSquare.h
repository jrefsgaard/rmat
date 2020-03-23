#ifndef HIST_LIKELIHOOD_CHI_SQUARE_H
#define HIST_LIKELIHOOD_CHI_SQUARE_H
#include <rmat/HistChiSquare.h>

namespace rmat {

class HistLikelihoodChiSquare : public HistChiSquare {
  public:
    HistLikelihoodChiSquare(std::shared_ptr<TH1> data, std::shared_ptr<TH1> model);
    ~HistLikelihoodChiSquare();
    
    double Evaluate();
};
}//namespace rmat
#endif//HIST_LIKELIHOOD_CHI_SQUARE_H
