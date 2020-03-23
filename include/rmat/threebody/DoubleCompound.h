#ifndef DOUBLE_COMPOUND_H
#define DOUBLE_COMPOUND_H
#include <vector>
#include <string>
#include <complex>
#include <armadillo>
#include "../Level.h"
#include "../Channel.h"
#include "../CompoundSystem.h"

namespace rmat {
namespace threebody {

class DoubleCompound {
  private:
    //Container for the primary and secondary compound system.
    std::vector<CompoundSystem> systems;
    
    //We must keep a list of standard 2b-channels for the primary decay.
    std::vector<Channel> channels;
    
    void CheckID(int sysid);
        
  public:
    DoubleCompound(int A1, int Z1, int Nl1, int A2, int Z2, int Nl2);
    ~DoubleCompound();
    
    arma::Mat<std::complex<double>> LevelMatrix(int sysid, double Ex, double J);
    
    void AddChannel(int sysid, Channel c, double threshold);
    void AddChannel(int sysid, std::string type, int l, double r0, double threshold);
    void ClearChannels();
    void UseInterpolation();

    void SetEnergy(int sysid, int lambda, double E);
    void SetJ(int sysid, int lambda, double J);
    void SetPartialWidths(int sysid, int lambda, std::vector<double> Gamma);
    void SetWidthAmplitude(int sysid, int lambda, int c, double gamma);
    void SetWidthAmplitudes(int sysid, int lambda, std::vector<double> gamma);
    void SetReducedWidth(int sysid, int lambda, int c, double gamma2);
    void SetReducedWidths(int sysid, int lambda, std::vector<double> gamma2);
    void SetDimensionlessWidth(int sysid, int lambda, int c, double theta2);
    void SetDimensionlessWidths(int sysid, int lambda, std::vector<double> theta2);
    
    double GetEnergy(int sysid, int lambda);
    double GetChannelEnergy(int sysid, double Ex, int c);
    double GetPartialWidth(int sysid, int lambda, int c);   //Observed partial width.
    double GetWidthAmplitude(int sysid, int lambda, int c); //Formal reduced width amplitude.
    double GetReducedWidth(int sysid, int lambda, int c);   //Formal reduced width.
    double GetDimensionlessWidth(int sysid, int lambda, int c);
    double GetWidth(int sysid, int lambda);                 //Total observed width.
    double GetRenormalisation(int sysid, int lambda); 
    int GetNLevels(int sysid);
    int GetNChannels(int sysid);
    std::vector<double> GetJs(int sysid);
    std::vector<int> GetLevelIndices(int sysid, double J);
    CompoundSystem & GetSystem(int sysid);

    //Caller must keep track of the channel ID's. The channels are ordered as
    //{(l,lambda),(l,lambda'),...,(l',lambda),(l',lambda'),...}.
    void LoadFcns(int channel, std::string datafile);
    //And a more user friendly version.
    void LoadFcns(int iL, int ilambda, std::string datafile);
    
    //Provides the channel ID based in L-index and lambda-index.
    int GetChannelIndex(int iL, int ilambda);
    
    //These are the unmodified primary channels.
    std::vector<Channel> & GetPrimaryChannels();
    
    //These are ordinary two-body channels.
    std::vector<Channel> & GetSecondaryChannels();
    
    int GetA(int sysid);
    int GetZ(int sysid);
};
} //namespace threebody;
} //namespace rmat;
#endif //DOUBLE_COMPOUND_H
