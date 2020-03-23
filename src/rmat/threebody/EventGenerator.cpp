#include <TMath.h>
#include <TGenPhaseSpace.h>
//#include <ausa/constants/Mass.h>
//#include <rmat/Nucleus.h>
#include <rmat/threebody/EventGenerator.h>

using namespace std;
using namespace TMath;

namespace rmat::threebody {

EventGenerator::EventGenerator() : He4{4,2,0,1}
{
  Q = 1000.;
}

EventGenerator::~EventGenerator() {}

void EventGenerator::SetQ(double q)
{
  Q = q;
}

double EventGenerator::Generate()
{
  return GenerateUniform();
}

double EventGenerator::GenerateUniform()
{
  //We'll need the mass of the alpha-particle.
  //Nucleus He4(4,2,0,1);

  double initialMass = (3 * He4.M() + Q)/1.e6; //Remember, GeV.
  TLorentzVector mother(0.0,0.0,0.0,initialMass);
  double masses[3] = {He4.M()/1e6, He4.M()/1e6, He4.M()/1e6}; //Remember, GeV.

  //Generate the decay.
  TGenPhaseSpace event;
  event.SetDecay(mother, 3, masses);
  double w = event.Generate();

  //Then we store the results in a safe place.
  products.at(0) = (*(event.GetDecay(0)))*1e6; //We like the results in keV.
  products.at(1) = (*(event.GetDecay(1)))*1e6;
  products.at(2) = (*(event.GetDecay(2)))*1e6;

  cmEnergies.at(0) = (products.at(0).Energy() - He4.M());
  cmEnergies.at(1) = (products.at(1).Energy() - He4.M());
  cmEnergies.at(2) = (products.at(2).Energy() - He4.M());

  return w;
}

const TLorentzVector & EventGenerator::GetProduct(int i)
{
  return products.at(i);
}


double EventGenerator::GetCmEnergy(int i)
{
  return cmEnergies.at(i);
}

array<TLorentzVector,3> EventGenerator::GetProducts()
{
  return products;
}

SimEvent EventGenerator::GetEvent()
{
  double x = Sqrt(3) * (GetCmEnergy(0) - GetCmEnergy(2)) / Q;
  double y = (2*GetCmEnergy(1) - GetCmEnergy(0) - GetCmEnergy(2))/Q;
  SimEvent event;
  event.decay = GetProducts();
  event.Q = Q;
  event.x = x;
  event.y = y;
  return event;
}

double EventGenerator::GetMaxWeight()
{
  return 0.5;
}
}//namespace rmat::threebody
