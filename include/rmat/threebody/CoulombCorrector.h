#ifndef COULOMB_CORRECTOR_H
#define COULOMB_CORRECTOR_H
#include <vector>
#include "rmat/Channel.h"
#include "DoubleCompound.h"

namespace rmat {
namespace threebody {

class CoulombCorrector {
  private:
    std::vector<Channel> pChannels;
    std::vector<Channel> sChannels;

  public:
    CoulombCorrector();
    CoulombCorrector(DoubleCompound &scheme, double radius);
    ~CoulombCorrector();

    double GetFSCI(double E, double E23, double E12, double E13, int pIndex, int sIndex);
};
} //namespace rmat
} //namespace threebody
#endif //COULOMB_CORRECTOR_H
