#include <rmat/threebody/DecayWeight.h>

using namespace std;

namespace rmat::threebody {

void DecayWeight::ExcludeQRange(double Qmin, double Qmax)
{
  array<double,2> range = {Qmin, Qmax};
  excludedQs.push_back(range);
}

vector<array<double,2>> & DecayWeight::GetExcludedQRanges()
{
  return excludedQs;
}

}
