#ifndef OBSERVED_SPECTRUM_H
#define OBSERVED_SPECTRUM_H
#include <memory>
#include <armadillo>
#include "Spectrum.h"
#include <TMatrixD.h>


namespace rmat {

class ObservedSpectrum : public Spectrum {
  /**
  * The ObservedSpectrum class allows the user to define a response matrix
  * for his experiment. When supplied with a model for the theoretical spectrum
  * this class will calculate the observed spectrum.
  */
  private:
    std::shared_ptr<Spectrum> model;
    arma::Mat<double> response;
    arma::Col<double> raw;
    arma::Col<double> observed;
    double emin, emax;
    int rebin;
    
    void SetResponse(TMatrixD &r); 
    void CalculateSpectrum();
    double Bin2Energy(int i);  //This is for calculating the raw spectrum
    int Energy2Bin(double E);  //This is for retrieving the correct element of the observed spectrum.
    void RebinSpectrum();
    
  public:
    ObservedSpectrum(std::shared_ptr<Spectrum> m, TMatrixD &r, double Emin, double Emax);
    ~ObservedSpectrum();
    
    double Strength(double Ec);
    void SetParameters(std::vector<double> par);
    void PrintParameters();
    int NDim();
   
    //Has only been tested with nice numbers that divide the response matrix exactly.
    void Rebin(int N = 2);
};
} //namespace rmat
#endif //OBSERVED_SPECTRUM_H
