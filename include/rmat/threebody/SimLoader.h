#ifndef SIM_LOADER_H
#define SIM_LOADER_H
#include <vector>
#include <string>
#include <TChain.h>
#include "rmat/threebody/SimEvent.h"

namespace rmat {
namespace threebody {

class SimLoader {
  private:
    std::string treename;
    //TODO: Could probably add support for alternative branch names.
    
  public:
    SimLoader(std::string treeName);
    ~SimLoader();
    
    /**
    * Extract a vector of SimEvents from an already loaded TChain.
    */
    std::vector<SimEvent> Load(TChain & data);
    
    /**
    * Extract a vector of SimEvents from the specified files. The fileName 
    * should contain the file path in the same format as that provided to
    * sprintf.
    */
    std::vector<SimEvent> Load(std::string filePath, int first = 0, int last = 0);
    
    void SetTreeName(std::string treeName);
    std::string GetTreeName();
};
} //namespace threebody
} //namespace rmat
#endif //SIM_LOADER_H
