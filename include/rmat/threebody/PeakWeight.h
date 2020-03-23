#ifndef PEAK_WEIGHT_H
#define PEAK_WEIGHT_H
#include "SimEvent.h"
#include "DoubleCompound.h"
#include "DecayWeight.h"
#include "ThreeBodyWeight.h"

namespace rmat {
namespace threebody {

class PeakWeight : public ThreeBodyWeight {
  private:
    int peakID;
    
  public:
    PeakWeight(std::string type, DoubleCompound s, int peak = 0);
    ~PeakWeight();
    
    //virtual double Calculate(std::array<TLorentzVector,3> &) override; 
    virtual double Calculate(SimEvent &) override;  
    virtual PeakWeight * Clone();

    void SetPeakID(int peak);
    int GetPeakID();
};
} //namespace rmat;
}//naespace threebody;
#endif //PEAK_WEIGHT_H
