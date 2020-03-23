#include <iostream>
#include <stdlib.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <rmat/threebody/Binner.h>

using namespace std;

namespace rmat::threebody {
    
Binner::Binner(string name, string title) : _name(name), _title(title)
{
  nX = 0;
  nY = 0;
  nZ = 0;
}

Binner::~Binner() {}
    
void Binner::SetX(int n, double min, double max)
{
  nX = n;
  minX = min;
  maxX = max;
}

void Binner::SetY(int n, double min, double max)
{
  nY = n;
  minY = min;
  maxY = max;
}

void Binner::SetZ(int n, double min, double max)
{
  nZ = n;
  minZ = min;
  maxZ = max;
}

int Binner::GetNX()
{
  return nX;
}

double Binner::GetMinX()
{
  return minX;
}

double Binner::GetMaxX()
{
  return maxX;
}

int Binner::GetNY()
{
  return nY;
}

double Binner::GetMinY()
{
  return minY;
}

double Binner::GetMaxY()
{
  return maxY;
}

int Binner::GetNZ()
{
  return nZ;
}

double Binner::GetMinZ()
{
  return minZ;
}

double Binner::GetMaxZ()
{
  return maxZ;
}
    
shared_ptr<TH1> Binner::Bin(vector<SimEvent> & events)
{
  shared_ptr<TH1> histogram;
  if(nX > 1 && nY > 1 && nZ > 1){
    //Standard 3D-spectrum with (x,y,Q)-coordinates.
    histogram = make_shared<TH3D> (_name.c_str(),_title.c_str(),nX,minX,maxX,nY,minY,maxY,nZ,minZ,maxZ);
    auto casted_histogram = dynamic_pointer_cast<TH3D>(histogram);
    for(SimEvent ev : events){ casted_histogram->Fill(ev.x,ev.y,ev.Q); }
  }
  else if(nX > 1 && nY > 1 && nZ <= 1){
    //Standard 2D-spectrum with (x,y)-values.
    histogram = make_shared<TH2D> (_name.c_str(),_title.c_str(),nX,minX,maxX,nY,minY,maxY);
    auto casted_histogram = dynamic_pointer_cast<TH2D>(histogram);
    for(SimEvent ev : events){ casted_histogram->Fill(ev.x,ev.y); }
  }
  else if(nX > 1 && nY <= 1 && nZ <= 1){
    //Standard 1D-spectrum with Q-values.
    histogram = make_shared<TH1D> (_name.c_str(),_title.c_str(),nX,minX,maxX);
    for(SimEvent ev : events){ histogram->Fill(ev.Q); }
  }
  else{
    cout << "  Binner::Bin(): Unknown binning requested:" ;
    cout << " nx = " << nX << ", ny = " << nY << ", nz = " << nZ << endl;
    exit(EXIT_FAILURE);
  }
  return histogram;  
}
    
shared_ptr<TH1> Binner::Bin(vector<SimEvent> & events, vector<bool> & accept)
{
  shared_ptr<TH1> histogram;
  if(events.size() != accept.size()){
    cout << "  Binner::Bin(): Event vector and acceptance vector have different sizes!" << endl;
    exit(EXIT_FAILURE);
  }
  int N = events.size();
  if(nX > 1 && nY > 1 && nZ > 1){
    //Standard 3D-spectrum with (x,y,Q)-coordinates.
    histogram = make_shared<TH3D> (_name.c_str(),_title.c_str(),nX,minX,maxX,nY,minY,maxY,nZ,minZ,maxZ);
    auto casted_histogram = dynamic_pointer_cast<TH3D>(histogram);
    for(int i=0; i<N; i++){
      if(accept[i]){
        SimEvent &ev = events[i];
        casted_histogram->Fill(ev.x,ev.y,ev.Q);
      }
    }
  }
  else if(nX > 1 && nY > 1 && nZ <= 1){
    //Standard 2D-spectrum with (x,y)-values.
    histogram = make_shared<TH2D> (_name.c_str(),_title.c_str(),nX,minX,maxX,nY,minY,maxY);
    auto casted_histogram = dynamic_pointer_cast<TH2D>(histogram);
    for(int i=0; i<N; i++){
      if(accept[i]){
        SimEvent &ev = events[i];
        //cout << "accepted, (x,y,) = (" << ev.x << "," << ev.y << ")" << endl;
        casted_histogram->Fill(ev.x,ev.y);
      }
    }
  }
  else if(nX > 1 && nY <= 1 && nZ <= 1){
    //Standard 1D-spectrum with Q-values.
    histogram = make_shared<TH1D> (_name.c_str(),_title.c_str(),nX,minX,maxX);
    for(int i=0; i<N; i++){
      if(accept[i]){
        SimEvent &ev = events[i];
        histogram->Fill(ev.Q);
      }
    }
  }
  else{
    cout << "  Binner::Bin(): Unknown binning requested:" ;
    cout << " nx = " << nX << ", ny = " << nY << ", nz = " << nZ << endl;
    exit(EXIT_FAILURE);
  }
  return histogram;
}
}//namespace rmat::threebody
