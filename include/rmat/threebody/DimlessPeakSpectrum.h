#ifndef DIMLESS_PEAK_SPECTRUM_H
#define DIMLESS_PEAK_SPECTRUM_H
#include <vector>
#include <string>
#include <rmat/NormalisedBetaSpectrum.h>
#include <rmat/CompoundSystem.h>
#include <rmat/Channel.h>

namespace rmat {
namespace threebody {
/**
 * A class for calculating the beta-delayed alpha spectrum for sequential decays
 * through a narrow resonance, for instance the 8Be(gs) peak. The spectrum is
 * normalised, so it uses the proper B(GT) values as parameters. The reduces
 * widths are dimensionless, to avoid issues with unphysical partial widths.
 */
class DimlessPeakSpectrum : public NormalisedBetaSpectrum {
  private:
    int peakID;
    CompoundSystem secondary_system;
    std::vector<Channel> channels_2b;
        
  public:
    /**
     * Constructor.
     * @param type Specifies whether it is a beta-plus or beta-minus decay.
     *             Possible values are "+" and "-".
     * @param primSys Primary compound system.
     * @param secSys Secondary compound system.
     */
    DimlessPeakSpectrum(std::string type, CompoundSystem primSys, CompoundSystem secSys);
    
    ~DimlessPeakSpectrum() = default;

    /**
     * Calculate the spectral strength/density.
     * @param Ec Channel energy, i.e. the energy above threshold.
     */
    virtual double Strength(double Ec);
    
    /** 
     * Calculate the spectral strength/density for a particular spin in the
     * primary compound system. Useful to visualise the contribution from
     * the different total spins, 0+, 1+, 2+, ...
     * @param Ec Channel energy, i.e. the energy above threshold.
     * @param J Spin of primary compound system
     */  
    virtual double Strength(double Ec, double J);
    
    /**
     * Set the parameters for the decay. The order of the parameters is
     *   {{E1,theta11,theta12,..,B(GT)1,E2,theta21,...},
     *   {E1,Gamma11,Gamma12,....,E2,Gamma21,...}}
     * where the first set of parameters is for the primary compound system
     * and the second set is for the secondary compound system.
     * Note: The theta's are actually theta^2 and the ordering is:
     * {(l,lambda),(l,lambda'),...,(l',lambda),(l',lambda'),...}.
     * @param par Vector containing the parameters.
     */
    virtual void SetParameters(std::vector<double> par);
    
    /**
     * Specify which peak in the secondary system constitutes the gate.
     * @param i Peak ID, i.e. 0 for the first level, 1 for the second and so on.
     */
    void SetPeakID(int i);
    int GetPeakID();
    
    CompoundSystem & GetSecondarySystem();
    
    std::vector<Channel> & GetTwoBodyChannels();
    
    virtual void PrintParameters();
    
    virtual int NDim();
    
};
}//threebody
}//rmat
#endif
