#ifndef EVENT_GENERATOR_H
#define EVENT_GENERATOR_H
#include <array>
#include <TLorentzVector.h>
#include <rmat/Nucleus.h>
#include <rmat/threebody/SimEvent.h>

namespace rmat::threebody {

class EventGenerator {

  private:
    double Q;
    std::array<TLorentzVector,3> products;  //Decay products
    std::array<double,3> cmEnergies;
    
    Nucleus He4;
    double GenerateUniform();

  public:
    /**
    * The EventGenerator will per default generate triple-alpha decays
    */
    EventGenerator();
    ~EventGenerator();

    /**
    * Set the available energy for the triple-alpha breakup. Units: keV.
    */
    void SetQ(double);

    /**
    * Generates a decay and returns the weight of that decay. The result of the decay
    * can be retrieved through the GetProduct() and GetCmEnergy() methods.
    */
    double Generate();
  
    std::array<TLorentzVector,3> GetProducts();
    const TLorentzVector & GetProduct(int);
    
    SimEvent GetEvent();

    /**
    * Returns the center of mass energies of the three decay products.
    */
    double GetCmEnergy(int);
    
    double GetMaxWeight();
};
}//namespace rmat::threebody
#endif
