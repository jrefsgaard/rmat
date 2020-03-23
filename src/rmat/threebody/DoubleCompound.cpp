#include <iostream>
#include <stdlib.h>  //exit(), EXIT_FAILURE
#include <rmat/threebody/DoubleCompound.h>

using namespace std;

namespace rmat::threebody {

DoubleCompound::DoubleCompound(int A1, int Z1, int Nl1, int A2, int Z2, int Nl2)
{
  CompoundSystem s1(A1,Z1,Nl1);
  CompoundSystem s2(A2,Z2,Nl2);
  systems.push_back(move(s1));
  systems.push_back(move(s2));
}

DoubleCompound::~DoubleCompound() {}

arma::Mat<std::complex<double>> DoubleCompound::LevelMatrix(int sysid, double Ex, double J)
{
  CheckID(sysid);
  return systems.at(sysid).LevelMatrix(Ex,J);
}

void DoubleCompound::AddChannel(int sysid, Channel c, double threshold)
{
  CheckID(sysid);
  
  if(sysid == 0){
    int Nb = systems.at(1).GetNLevels();
    for(int i=0; i<Nb; i++){
      Channel ci = c;  //Make Nb copies
      systems.at(0).AddChannel(move(ci),threshold);
    }
    Channel cs = c;
    channels.push_back(cs);
  }
  else{
    systems.at(1).AddChannel(move(c),threshold);
  }
}

void DoubleCompound::AddChannel(int sysid, std::string type, int l, double r0, double threshold)
{
  CheckID(sysid);
  
  if(sysid == 0){
    int Nb = systems.at(1).GetNLevels();
    for(int i=0; i<Nb; i++){
      systems.at(0).AddChannel(type,l,r0,threshold);
    }
    Channel cs = systems.at(0).GetChannels().back();
    channels.push_back(cs);
  }
  else{
    systems.at(1).AddChannel(type,l,r0,threshold);
  }
}

void DoubleCompound::ClearChannels()
{
  systems.at(0).ClearChannels();
  systems.at(1).ClearChannels();
  channels.clear();
}

void DoubleCompound::UseInterpolation()
{
  systems.at(0).UseInterpolation();
  systems.at(1).UseInterpolation();
  for(auto & c : channels) c.UseInterpolation();
}

void DoubleCompound::SetWidthAmplitude(int sysid, int lambda, int c, double gamma)
{
  CheckID(sysid);
  systems.at(sysid).SetWidthAmplitude(lambda,c,gamma);
}

void DoubleCompound::SetWidthAmplitudes(int sysid, int lambda, std::vector<double> gamma)
{
  CheckID(sysid);
  systems.at(sysid).SetWidthAmplitudes(lambda,gamma);
}

void DoubleCompound::SetReducedWidth(int sysid, int lambda, int c, double gamma2)
{
  CheckID(sysid);
  systems.at(sysid).SetReducedWidth(lambda,c,gamma2);
}

void DoubleCompound::SetReducedWidths(int sysid, int lambda, std::vector<double> gamma2)
{
  CheckID(sysid);
  systems.at(sysid).SetReducedWidths(lambda,gamma2);
}

void DoubleCompound::SetDimensionlessWidth(int sysid, int lambda, int c, double theta2)
{
  CheckID(sysid);
  systems.at(sysid).SetDimensionlessWidth(lambda,c,theta2);
} 

void DoubleCompound::SetDimensionlessWidths(int sysid, int lambda, vector<double> theta2)
{
  CheckID(sysid);
  systems.at(sysid).SetDimensionlessWidths(lambda,theta2);
}

void DoubleCompound::SetPartialWidths(int sysid, int lambda, std::vector<double> Gamma)
{
  CheckID(sysid);
  systems.at(sysid).SetPartialWidths(lambda,Gamma);
}

void DoubleCompound::SetEnergy(int sysid, int lambda, double E)
{
  CheckID(sysid);
  systems.at(sysid).SetEnergy(lambda,E);
}

void DoubleCompound::SetJ(int sysid, int lambda, double J)
{
  CheckID(sysid);
  systems.at(sysid).GetLevel(lambda).SetJ(J);
}

void DoubleCompound::CheckID(int sysid)
{
  if(!(sysid == 0 || sysid == 1)){
    cout << "  DoubleCompound::CheckID(): Wrong sysid: "<< sysid << "!" << endl;
    exit(EXIT_FAILURE);
  }
}

double DoubleCompound::GetEnergy(int sysid, int lambda)
{
  CheckID(sysid);
  return systems.at(sysid).GetEnergy(lambda);
}

double DoubleCompound::GetChannelEnergy(int sysid, double Ex, int c)
{
  CheckID(sysid);
  return systems.at(sysid).GetChannelEnergy(Ex,c);
}

double DoubleCompound::GetPartialWidth(int sysid, int lambda, int c)   //Observed partial width.
{
  CheckID(sysid);
  return systems.at(sysid).GetPartialWidth(lambda,c);
}

double DoubleCompound::GetWidthAmplitude(int sysid, int lambda, int c) //Formal reduced width amplitude.
{
  CheckID(sysid);
  return systems.at(sysid).GetWidthAmplitude(lambda,c);
}

double DoubleCompound::GetReducedWidth(int sysid, int lambda, int c)   //Formal reduced width.
{
  CheckID(sysid);
  return systems.at(sysid).GetReducedWidth(lambda,c);
}

double DoubleCompound::GetDimensionlessWidth(int sysid, int lambda, int c)
{
  CheckID(sysid);
  return systems.at(sysid).GetDimensionlessWidth(lambda,c);
}

double DoubleCompound::GetWidth(int sysid, int lambda)                 //Total observed width.
{
  CheckID(sysid);
  return systems.at(sysid).GetWidth(lambda);
}

double DoubleCompound::GetRenormalisation(int sysid, int lambda)
{
  CheckID(sysid);
  return systems.at(sysid).GetRenormalisation(lambda);
} 

int DoubleCompound::GetNLevels(int sysid)
{
  CheckID(sysid);
  return systems.at(sysid).GetNLevels();
}

int DoubleCompound::GetNChannels(int sysid)
{
  CheckID(sysid);
  return systems.at(sysid).GetNChannels();
}

vector<double> DoubleCompound::GetJs(int sysid)
{
  CheckID(sysid);
  return systems.at(sysid).GetJs();
}

vector<int> DoubleCompound::GetLevelIndices(int sysid, double J)
{
  CheckID(sysid);
  return systems.at(sysid).GetLevelIndices(J);
}

void DoubleCompound::LoadFcns(int channel, string datafile)
{
  systems.at(0).GetChannel(channel).UseInterpolation(datafile);
}

void DoubleCompound::LoadFcns(int iL, int ilambda, std::string datafile)
{
  int index = GetChannelIndex(iL,ilambda);
  LoadFcns(index,datafile);
}

int DoubleCompound::GetChannelIndex(int iL, int ilambda)
{
  int NL = channels.size();
  int Nb = systems.at(1).GetNLevels();
  return iL * Nb + ilambda;
}

vector<Channel> & DoubleCompound::GetPrimaryChannels()
{
  return channels;
}

vector<Channel> & DoubleCompound::GetSecondaryChannels()
{
  return systems.at(1).GetChannels();
}

int DoubleCompound::GetA(int sysid)
{
  CheckID(sysid);
  return systems.at(sysid).A();
}

int DoubleCompound::GetZ(int sysid)
{
  CheckID(sysid);
  return systems.at(sysid).Z();
}

CompoundSystem & DoubleCompound::GetSystem(int sysid)
{
  CheckID(sysid);
  return systems.at(sysid);
}
}//namespace rmat::threebody
