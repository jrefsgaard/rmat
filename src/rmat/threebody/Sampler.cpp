#include <iostream>
#include <stdlib.h>
#include <rmat/threebody/Sampler.h>
    
using namespace std;

namespace rmat::threebody {

Sampler::Sampler() : rGen{0} {}

Sampler::~Sampler() {}

vector<SimEvent> Sampler::Sample(vector<SimEvent> & raw, vector<double> & weights, double maxWeight)
{
  if(raw.size() != weights.size()){
    cout << "  Sampler::Sample(): Event vector and weight vector are of different sizes!" << endl;
    exit(EXIT_FAILURE);
  }
  int N = raw.size();
  vector<SimEvent> result;
  for(int i=0; i<N; i++){
    if(Accept(weights.at(i),maxWeight)) result.push_back(raw.at(i));
  }
  
  return result;
}

vector<bool> Sampler::Sample(vector<double> & weights, double maxWeight)
{
  int N = weights.size();
  vector<bool> result (N,false);
  for(int i=0; i<N; i++){
    result[i] = Accept(weights.at(i),maxWeight);
  }
  
  return result;   
}

long long int Sampler::Sample(vector<double> & weights, double maxWeight, vector<bool> & accepted)
{
  if(weights.size() != accepted.size()){
    cout << "  Sampler::Sample(): Weight vector and accept vector are of different sizes!" << endl;
    exit(EXIT_FAILURE);
  }
  long long int N = weights.size();
  long long int Naccept = 0;
  for(int i=0; i<N; i++){
    bool ai = Accept(weights.at(i),maxWeight);
    accepted[i] = ai;
    if(ai) Naccept++;
  }
  return Naccept;
}

bool Sampler::Accept(double weight, double maxWeight)
{
  if(rGen.Uniform(0,maxWeight) < weight) return true; //Event accepted.
  else return false; //Event rejected.
}
}//namespace rmat::threebody
