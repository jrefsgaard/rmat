#ifndef LEVEL_H
#define LEVEL_H
#include <vector>

namespace rmat {

class Level {
  /*
  * The Level class stores basic information on a nuclear level/state/resonance.
  * The energy should be given as energy above the ground state.
  */
  
  private:
    double energy;
    std::vector<double> widths;
    int jx2;
    int pi;
    
  public:
    Level(double E, std::vector<double> gamma, double j = 0., int parity = 1);
    ~Level();

    void SetE(double E);
    void SetEnergy(double E);
    void SetJ(double);
    void SetParity(int);
    void SetWidthAmplitude(int c, double gammac);
    void SetGammas(std::vector<double> gammac);

    double E() const;
    double GetEnergy() const;
    double J() const;
    int Parity() const;
    int Pi() const;
    std::vector<double> & Gammas();
    std::vector<double> & GetWidthAmplitudes();
    double GetWidthAmplitude(int c);
};
} //namespace rmat;
#endif  //LEVEL_H
