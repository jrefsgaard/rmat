#include <iostream>
#include <TDatime.h>
#include <rmat/threebody/Simulator.h>
#include <rmat/threebody/WeightCalculator.h>
#include <rmat/threebody/Sampler.h>
#include <rmat/threebody/Binner.h>

using namespace std;

namespace rmat::threebody {

Simulator::Simulator(shared_ptr<DecayWeight> w, vector<SimEvent> pool)
: weight(w), eventPool(move(pool)), accepted(eventPool.size(),false)
{}
    
Simulator::~Simulator() {}
    
long long int Simulator::Simulate()
{
  TDatime timer;
  
  WeightCalculator calculator(weight);
  //cout << "Simulator: Starting weight calculation: " << endl; timer.Set(); timer.Print();
  vector<double> weights = calculator.CalculateWeights(eventPool);
  //cout << "Simulator: Calculation done: " << endl; timer.Set();timer.Print();
  
  double maxWeight = calculator.MaxWeight();
  
  //cout << weights.size() << " weights calculated, max weight = " << maxWeight << endl;
  
  //cout << "Simulator: Sampling started: " << endl; timer.Set(); timer.Print();
  Sampler sampler;
  long long int N = sampler.Sample(weights,maxWeight,accepted); 
  //cout << "Simulator: Sampling finished: " << endl; timer.Set(); timer.Print(); 
  
  return N;
}


shared_ptr<TH1> Simulator::GetSpectrum(int nx, double xMin, double xMax,
                                       int ny, double yMin, double yMax,
                                       int nz, double zMin, double zMax)
{
  //cout << "line 44" << endl;
  Binner binner("simulation","Result from 3body-simulation");
  //cout << "line 46" << endl;
  binner.SetX(nx,xMin,xMax);
  binner.SetY(ny,yMin,yMax);
  binner.SetZ(nz,zMin,zMax);
  
  return binner.Bin(eventPool,accepted);
}

void Simulator::SetPool(std::vector<SimEvent> pool)
{
  eventPool = pool;
}
 
void Simulator::SetDecayWeight(shared_ptr<DecayWeight> w)
{
  weight = w;
}

vector<SimEvent> & Simulator::GetPool()
{
  return eventPool;
}

DecayWeight & Simulator::GetDecayWeight()
{
  return *(weight.get());
}

unsigned int Simulator::NDim()
{
  return weight->NDim();
}

void Simulator::SetParameters(vector<double> par)
{
  weight->SetParameters(par);
}
} //namespace rmat::threebody
