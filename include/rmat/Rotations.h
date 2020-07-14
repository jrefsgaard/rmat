#ifndef ROTATIONS_H
#define ROTATIONS_H
#include <TRotation.h>
#include <vector>

/**
 * Class providing a table of rotations randomly distributed on the unit sphere.
 */
class Rotations {
  private:
    int n_rot;
    std::vector<TRotation> rotations;
    
  public:
    Rotations(int N = 100);
    ~Rotations() = default;

    const TRotation & GetRotation(int i);
    
    std::vector<TRotation> & GetRotations();

    //std::vector<TRotation> GetRotation(int i);
};
#endif
