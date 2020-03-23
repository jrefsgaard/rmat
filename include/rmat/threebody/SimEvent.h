#ifndef SIM_EVENT_H
#define SIM_EVENT_H
#include <array>
#include <TLorentzVector.h>

namespace rmat::threebody {

struct SimEvent {
  std::array<TLorentzVector,3> decay;
  double Q;
  double x;
  double y;
};
}//namespace rmat::threebody;
#endif //SIM_EVENT_H
