#include <iostream>
#include <rmat/threebody/CheatSimulator.h>
#include <rmat/threebody/DalitzPlot.h>

using namespace std;

namespace rmat::threebody {

CheatSimulator::CheatSimulator(shared_ptr<DecayWeight> w,
                               shared_ptr<TH3D> simulated,
                               shared_ptr<TH3D> accepted )
: Simulator(w)
{
  shared_ptr<DalitzPlot> dplot (new DalitzPlot(w));
  plot = shared_ptr<ObservedDalitzPlot>
               (new ObservedDalitzPlot(dplot,simulated,accepted));
}

 
double CheatSimulator::Simulate()
{
  result = plot->GetSpectrum();
  double N = result->Integral();
  return N;
}


shared_ptr<TH1> CheatSimulator::GetSpectrum(int nx, double xMin, double xMax,
                                       int ny, double yMin, double yMax,
                                       int nz, double zMin, double zMax)
{
  return result;
}

} //namespace rmat::threebody
