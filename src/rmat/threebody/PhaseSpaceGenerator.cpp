#include <rmat/threebody/PhaseSpaceGenerator.h>

using namespace std;

namespace rmat::threebody {

PhaseSpaceGenerator::PhaseSpaceGenerator() : rGen{0} {}

PhaseSpaceGenerator::~PhaseSpaceGenerator() {}
    
vector<SimEvent> PhaseSpaceGenerator::Generate(int N, double Q)
{
  vector<SimEvent> result;
  generator.SetQ(Q);
  double maxWeight = generator.GetMaxWeight();
  while(result.size() < N){
    double wi = generator.Generate();
    bool accept = sampler.Accept(wi,maxWeight);
    if(accept) result.push_back(generator.GetEvent());
  }
  
  return result;    
}
    
vector<SimEvent> PhaseSpaceGenerator::Generate(int N, double Qmin, double Qmax)
{
  vector<SimEvent> result;
  double maxWeight = generator.GetMaxWeight();
  while(result.size() < N){
    double qi = rGen.Uniform(Qmin,Qmax);
    generator.SetQ(qi);
    double wi = generator.Generate();
    bool accept = sampler.Accept(wi,maxWeight);
    if(accept) result.push_back(generator.GetEvent());
  }
  
  return result;  
}
}
