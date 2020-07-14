#include <rmat/Rotations.h>
#include <cstdlib>
#include <iostream>
#include <TVector3.h>
#include <TRandom3.h>
#include <TMath.h>

using namespace std;

Rotations::Rotations(int N) : n_rot(N)
{
  if(n_rot <= 0){
    cout << "Rotations::Rotations(): Cannot initialise less than one rotation." << endl;
    exit(EXIT_FAILURE);
  }
  
  TRandom3 rGen;
  for(int i=0; i<n_rot; i++){
    double x, y, z;
    rGen.Sphere(x,y,z,1.);
    TVector3 v(x,y,z);
    double t = rGen.Uniform(0.,2*TMath::Pi());
    TRotation r;
    r.Rotate(t,v);
    rotations.push_back(r);
  }
}

const TRotation & Rotations::GetRotation(int i)
{
  return rotations.at(i%n_rot);
}

vector<TRotation> & Rotations::GetRotations()
{
  return rotations;
}

/*
vector<TRotation> Rotations::GetRotation(int N)
{
  auto first = rotations.begin();
  auto last = rotations.begin() + N;
  vector<TRotation> vec(first,last);
  return vec;
}
*/
