#ifndef CHANNEL_H
#define CHANNEL_H
#include <string>
#include "Math/Interpolator.h"
#include "Interpolator.h"
#include "Pair.h"
#include "Constants.h"
//#include <memory>

namespace rmat {

class Channel {
  /**
  * This class will calculate some properties often used in R-matrix analysis.
  * It needs to know the constituent particles of the breakup,
  * relative angular momentum and the desired channel radius.
  * Units: Mass and energy in keV, lengths in fm.
  */

  private:
    Pair pair;       //Which particle pair does the breakup proceed through.
    int l;           //Orbital angular momentum between the fragments.
    double r;        //Channel radius in fm.

    //The following is for fast evaluation.
    bool interpolate;  //Do interpolation or not?
    double emax, emin; //Upper and lower energy limit of the existing interpolation. emin = -0.5*emax.
    double estep;      //Grid of 'exact' points on which to perform the interpolation.

    Interpolator Pl;      //Penetrability.
    Interpolator S;       //Shift function.
    Interpolator dS;      //Derivative of shift function.
    Interpolator Phi;     //Hard-sphere phase shift.

    void MakeInterpolation(double, double);  //Initialises the interpolation.
    void CheckAndExpand(double);

    double ExactPenetrability(double);    //Evaluation using high-precision Coulomb functions.
    double FastPenetrability(double);     //Evaluation using interpolation.

    double ExactShiftFunction(double);
    double FastShiftFunction(double);

    double ExactShiftDeriv(double);
    double FastShiftDeriv(double);

    double ExactHardSphere(double);
    double FastHardSphere(double);

  public:
    
    Channel(Pair pc = Pair(), int lc = 0, double rc = 1.42, bool interpc = false);
    ~Channel();

    void SetPair(Pair&);

    /**
    * Change the relative angular momentum between the breakup-fragments.
    */
    void SetL(int);
 
    /**
    * The channel radius may be varied. The default value is calculated using
    * r0 = 1.42fm, but here you can set the total radius.
    */
    void SetRadius(double);
    
    /**
    * Set the channel radius by specifying the nucleon radius (default r0 = 1.42fm).
    */
    void SetNucleonRadius(double r0);

    const Pair & GetPair();
    int L();
    double Radius();

    /**
    * The R-matrix functions may be evaluated faster using interpolation. In this
    * implementation we use cubic spline interpolation with user defined limits
    * and step size. If the maximum energy is later exceeded, the interpolation
    * is automatically expanded to incorporate a larger energy region. One should
    * be careful with the fast option, since it takes time to initially set up the
    * interpolation.
    */
    void UseInterpolation(bool interp = true, double Emax = 10000., double Estep = 20.);
    
    /**
    * We also allow the possibility of pre-calculated R-matrix functions to be 
    * used. The caller must supply a five-column txt-file, where the first
    * column is the channel energy (in keV) and the remaining columns contain
    * P, S, dS and phi. Be careful with the argument (compiler converts string 
    * literals to bool).
    */
    void UseInterpolation(std::string datafile);

    /**
    * The following are some common functions needed in R-matrix analysis. 
    * Arguments are all energy above the channel threshold in keV.
    */
    double Penetrability(double);
    double ShiftFunction(double);
    double ShiftFunctionDeriv(double);
    double CoulombShift(double);
    double HardSphereShift(double);
    double BarrierHeight();        //Height of Coulomb barrier in [keV].
    double Rho(double);            //'Distance' parameter, (wave number) * r.
    double WignerLimit();
};
}  //namespace rmat
#endif  //CHANNEL_H
