#ifndef CHI_SQUARE_H
#define CHI_SQUARE_H
#include <memory>
#include <vector>
#include <TH1.h>
#include <Math/IFunction.h>
#include "Spectrum.h"

namespace rmat {

class ChiSquare : public ROOT::Math::IMultiGenFunction {
  /**
  * ChiSquare is a type of error function which will compare a theoretical
  * model with an experimental data set. Since we inherit from IMultiGenFunction
  * it is possible to pass an object of the ChiSquare-class to ROOT's
  * Minimizer class.
  */
  private:
    //XXX
    virtual double DoEval(const double *x) const override;

  protected:
    std::shared_ptr<TH1> data;
    std::shared_ptr<Spectrum> model;
    mutable int nBins;
    double binning;
    double minBin, maxBin;
    
  public:
    ChiSquare(std::shared_ptr<TH1> d, std::shared_ptr<Spectrum> m);
    ~ChiSquare();
    
    //XXX
    virtual ROOT::Math::IBaseFunctionMultiDim *Clone() const override;
    //XXX
    virtual unsigned int NDim() const override;

    virtual double Evaluate() const;
    //It is probably bad practice, but I want to be able to evaluate the ChiSquare
    //also in a public method.
    void SetParameters(const double *par) const;
    
    void SetData(std::shared_ptr<TH1> d);
    void SetModel(std::shared_ptr<Spectrum> mo);
    void SetBinning(double b);
    void SetBinRange(int min, int max);
    void SetRange(double min, double max);
        
    TH1 & GetData();
    Spectrum & GetModel();
    double GetBinning();
    int GetMinBin();
    int GetMaxBin();
    double GetMin();
    double GetMax();
    int GetCounts();
    int Ndf();
    int GetNBins();
};
} //namespace rmat
#endif //CHI_SQUARE_H
