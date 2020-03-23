#ifndef HIST_SYMMETRIC_CHI_SQUARE_H
#define HIST_SYMMETRIC_CHI_SQUARE_H
#include <rmat/HistChiSquare.h>

namespace rmat {

class HistSymmetricChiSquare : public HistChiSquare {
  public:
    HistSymmetricChiSquare(std::shared_ptr<TH1> data, std::shared_ptr<TH1> model);
    ~HistSymmetricChiSquare();
    
    double Evaluate();
};
}//namespace rmat
#endif//HIST_SYMMETRIC_CHI_SQUARE_H
