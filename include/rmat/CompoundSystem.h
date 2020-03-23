#ifndef COMPOUND_SYSTEM_H
#define COMPOUND_SYSTEM_H
#include <vector>
#include <complex>
//#include <memory>
#include <string>
#include <armadillo>
#include "Channel.h"
#include "Level.h"

namespace rmat {

class CompoundSystem {
  private :
    int a;
    int z;
    std::vector<Level> levels;
    std::vector<Channel> channels;
    std::vector<double> thresholds;
    
    void CheckAndResize();
    void EnsureWidthSize(std::vector<double> &);
   
  public :
    CompoundSystem(int A, int Z, int Nl);
    ~CompoundSystem();
    
    arma::Mat<std::complex<double>> LevelMatrix(double Ex, double J);
    
    Channel & AddChannel(Channel c, double threshold);
    Channel & AddChannel(std::string type, int l, double r0, double threshold);
    void AddLevel(Level level);
    void AddLevel(double E, std::vector<double> gamma);
    void AddLevel(double E);
    void ClearChannels();
    void ClearLevels();
    void UseInterpolation(bool interp = true, double Emax = 10000., double Estep = 20.);

    void SetEnergy(int lambda, double E);
    void SetPartialWidths(int lambda, std::vector<double> Gamma);
    void SetWidthAmplitude(int lambda, int c, double gamma);
    void SetWidthAmplitudes(int lambda, std::vector<double> gamma);
    void SetReducedWidth(int lambda, int c, double gamma2);
    void SetReducedWidths(int lambda, std::vector<double> gamma2);
    void SetDimensionlessWidth(int lambda, int c, double theta2);
    void SetDimensionlessWidths(int lambda, std::vector<double> theta2);    
        
    double GetEnergy(int lambda);
    double GetChannelEnergy(double Ex, int c);
    double GetPartialWidth(int lambda, int c);   //Observed partial width.
    double GetWidthAmplitude(int lambda, int c); //Formal reduced width amplitude.
    double GetReducedWidth(int lambda, int c);   //Formal reduced width.
    double GetDimensionlessWidth(int lambda, int c);
    double GetWidth(int lambda);                 //Total observed width.
    double GetRenormalisation(int lambda);    
    Level & GetLevel(int lambda);
    std::vector<Level> & GetLevels();
    Channel & GetChannel(int c);
    std::vector<Channel> & GetChannels();
    double GetThreshold(int c);
    std::vector<double> & GetThresholds();
    int GetNLevels();
    int GetNChannels();
    int A();
    int Z();
    std::vector<double> GetJs();
    std::vector<int> GetLevelIndices(double J);
   
    const void PrintParameters();
};
} //namespace rmat
#endif //COMPOUND_SYSTEM_H
