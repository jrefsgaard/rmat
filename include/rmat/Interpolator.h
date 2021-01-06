#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H
#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include "Math/Interpolator.h"
#include "Math/InterpolationTypes.h"

namespace rmat {

class Interpolator {
  /**
  * This class provides a copyable version of ROOT's otherwise
  * excellent Interpolator class, which has private copy- and
  * copy-assignment constructors. Since the deep copy takes time,
  * however, it should not be copied around needlessly.
  */
  
  private:
    std::vector<double> x, y;
    ROOT::Math::Interpolation::Type type;
    std::unique_ptr<ROOT::Math::Interpolator> interpolator;

    //ROOT::Math::Interpolator does not handle zero-length input very well,
    //so we take special care.    
    std::unique_ptr<ROOT::Math::Interpolator> SafeConstructor
      (const std::vector<double>& vx, const std::vector<double>& vy, ROOT::Math::Interpolation::Type t = ROOT::Math::Interpolation::kCSPLINE);
      
    bool CheckBounds(double var) const;
    
  public:
    Interpolator(unsigned int ndata = 0, ROOT::Math::Interpolation::Type t = ROOT::Math::Interpolation::kCSPLINE);
    Interpolator(const std::vector<double>& vx, const std::vector<double>& vy, ROOT::Math::Interpolation::Type t = ROOT::Math::Interpolation::kCSPLINE);
    bool SetData(const std::vector<double>& vx, const std::vector<double>& vy);
    bool SetData(unsigned int ndata, const double* vx, const double* vy);
    
    //Copy constructor.
    Interpolator(const Interpolator &);
    //Copy-assignment operator.
    Interpolator & operator = (const Interpolator &);
    //Move constructor.
    Interpolator(Interpolator&&);
    //Move-assignment operator
    Interpolator & operator = (Interpolator&&);
    
    //Reimplement all const-functions from ROOT::Math::Interpolator
    double Deriv(double var) const;// {return interpolator->Deriv(var);};
    double Deriv2(double var) const;// {return interpolator->Deriv2(var);};
    double Eval(double var) const;// {
      //if(var < x.front() || var > x.back())
      //  std::cout << "Evaluating at " << var << std::endl;
      //return interpolator->Eval(var);};
    double Integ(double a, double b) const;// {return interpolator->Integ(a,b);};
    std::string Type() const {return interpolator->Type();};
    std::string TypeGet() const {return interpolator->TypeGet();};
};
} //namespace rmat.
#endif //INTERPOLATOR_H
