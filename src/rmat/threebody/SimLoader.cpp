#include <iostream>
#include <TClonesArray.h>
#include <TBranch.h>
#include <TLorentzVector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <rmat/threebody/SimLoader.h>

using namespace std;

namespace rmat::threebody {
    
SimLoader::SimLoader(string treeName) : treename(treeName) {}

SimLoader::~SimLoader() {}
    
vector<SimEvent> SimLoader::Load(TChain & data)
{
  TTreeReader reader(&data);
  TTreeReaderValue<double> x(reader,"x");
  TTreeReaderValue<double> y(reader,"y");
  TTreeReaderValue<double> Q(reader,"Q");
  TTreeReaderArray<TLorentzVector> pAlpha(reader,"pAlpha");
  
  vector<SimEvent> events;
  while(reader.Next()){
    SimEvent ev;
    ev.decay[0] = pAlpha[0];
    ev.decay[1] = pAlpha[1];
    ev.decay[2] = pAlpha[2];
    ev.Q = *Q;
    ev.x = *x;
    ev.y = *y;
    events.push_back(ev);
  }
  
  return events;
}
    
vector<SimEvent> SimLoader::Load(string filePath, int first, int last)
{
  TChain chain(treename.c_str());
  
  //We add the files to the TChain.
  char fileName[256];
  for(int i=first; i<=last; i++){
    sprintf(fileName,filePath.c_str(),i);
    chain.Add(fileName);
  }
  
  return Load(chain);
}
    
void SimLoader::SetTreeName(string treeName)
{
  treename = treeName;
}

string SimLoader::GetTreeName()
{
  return treename;
}
} //namespace rmat::threebody
