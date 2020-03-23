#ifndef BINNER_H
#define BINNER_H
#include <vector>
#include <memory>
#include <string>
#include <TH1.h>
#include <rmat/threebody/SimEvent.h>

namespace rmat {
namespace threebody {

class Binner {
  private:
    std::string _name;
    std::string _title;
    int nX;
    double minX, maxX;
    int nY;
    double minY, maxY;
    int nZ;
    double minZ, maxZ;
    
  public:
    Binner(std::string name, std::string title);
    ~Binner();
    
    void SetX(int n, double min, double max);
    void SetY(int n, double min, double max);
    void SetZ(int n, double min, double max);
    int GetNX();
    double GetMinX();
    double GetMaxX();
    int GetNY();
    double GetMinY();
    double GetMaxY();
    int GetNZ();
    double GetMinZ();
    double GetMaxZ();
    
    std::shared_ptr<TH1> Bin(std::vector<SimEvent> & events);
    
    std::shared_ptr<TH1> Bin(std::vector<SimEvent> & events, std::vector<bool> & accept);
};
}//namespace threebody
}//namespace rmat
#endif //BINNER_H
