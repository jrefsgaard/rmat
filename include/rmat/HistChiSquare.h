#ifndef HIST_CHI_SQUARE_H
#define HIST_CHI_SQUARE_H
#include <memory>
#include <TH1.h>

namespace rmat {

class HistChiSquare {
  protected:
    std::shared_ptr<TH1> _data;
    std::shared_ptr<TH1> _model;
    int nBins;
    
    void CheckDimension();
    
  public:
    /**
    * Calculates the distance between two histograms. Since all ROOT-histograms
    * inherit from TH1, it will take all types and dimensions.
    */    
    HistChiSquare(std::shared_ptr<TH1> data, std::shared_ptr<TH1> model);
    ~HistChiSquare();
    
    void SetData(std::shared_ptr<TH1> data);
    void SetModel(std::shared_ptr<TH1> model);
    TH1 & GetData();
    TH1 & GetModel();
    
    virtual double Evaluate();
    int GetNBins();
};
}//namespace rmat
#endif //HIST_CHI_SQUARE_H
