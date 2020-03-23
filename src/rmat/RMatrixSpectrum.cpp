#include <rmat/RMatrixSpectrum.h>

using namespace std;

namespace rmat {

RMatrixSpectrum::RMatrixSpectrum(CompoundSystem s) : system(move(s)) {}

RMatrixSpectrum::~RMatrixSpectrum() {}

void RMatrixSpectrum::SetOutChannels(std::vector<int> channels)
{
  int Nc = system.GetNChannels();
  for(int ci : channels){
    if(ci < 0 || ci >= Nc){
      cout << "  RMatrixSpectrum::SetOutChannels(): Channel " << ci << " does not exist!" << endl;
      exit(EXIT_FAILURE);
    }
  }
  outChannels = channels;
}

void RMatrixSpectrum::SetOutChannel(int c)
{
  vector<int> cv = {c};
  SetOutChannels(cv);
}

const vector<int> & RMatrixSpectrum::GetOutChannels()
{
  return outChannels;
}

CompoundSystem & RMatrixSpectrum::GetSystem()
{
  return system;
}
}//namespace rmat
