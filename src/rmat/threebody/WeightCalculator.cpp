#include <iostream>
#include <omp.h>
#include <rmat/threebody/WeightCalculator.h>
#include <TDatime.h>

using namespace std;

namespace rmat::threebody {

WeightCalculator::WeightCalculator(shared_ptr<DecayWeight> w) : weight(w), maxWeight(0.) {}

WeightCalculator::~WeightCalculator() {}
    
vector<double> WeightCalculator::CalculateWeights(vector<SimEvent> & raw)
{
  maxWeight = 0.;
  
  int nDecays = raw.size();
  vector<double> weights (nDecays);
  vector<double> maxWeights;
  vector<unique_ptr<DecayWeight>> calculators;
  #pragma omp parallel// num_threads(8)
  {
    int nThreads = omp_get_num_threads();
    //We establish number of threads and clone the weight calculators.
    #pragma omp master
    {
      //cout << "line 27 " << nThreads << endl;
      maxWeights.resize(nThreads,-1);
      for(int i=0; i<nThreads; i++){
        unique_ptr<DecayWeight> ci = unique_ptr<DecayWeight>(weight->Clone());
        calculators.push_back(move(ci));
      }
      //TDatime timer;
      //cout << "Starting loop: "<< endl; timer.Set(); timer.Print();
    }
    //cout << "line 35,  calculators.size() = " << calculators.size() << endl;
    #pragma omp barrier
    #pragma omp for
    for(int i=0; i<nDecays; i++){
      int id = omp_get_thread_num();
      //cout << "line 40 " << id << "  calculators.size() = " << calculators.size() << endl;
      //double Wi = calculators.at(id)->Calculate(raw.at(i).decay);
      double Wi = calculators.at(id)->Calculate(raw.at(i));
      //cout << "Wi = " << Wi << endl;
      weights.at(i) = Wi;
      if(Wi > maxWeights.at(id)) maxWeights.at(id) = Wi;
    }
    #pragma omp barrier
    #pragma omp master
    {
      //TDatime timer;
      //cout << "Finishing loop: "<< endl; timer.Set(); timer.Print(); 
      for(int i=0; i<nThreads; i++){
        if(maxWeights.at(i) > maxWeight){ maxWeight = maxWeights.at(i);}
      }
    }
  }
  
  return weights;
}
    
double WeightCalculator::MaxWeight()
{
  return maxWeight;
}

void WeightCalculator::SetWeight(std::shared_ptr<DecayWeight> w)
{
  weight = w;
}

DecayWeight & WeightCalculator::GetWeight()
{
  return *(weight.get());
}
    
}//namespace rmat::threebody;
