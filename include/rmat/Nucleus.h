#ifndef NUCLEUS_H
#define NUCLEUS_H
#include "Level.h"

namespace rmat {

class Nucleus {
  /**
  * The Nucleus object is a container to hold basic information describing
  * a nucleus.
  */

  private:
    int a;
    int q;
    double m;
    int jx2;
    int pi;

  public:
    /**
    * The default constructor will produce a 4He-nucleus.
    */
    Nucleus(int A = 4, int Z = 2, double J = 0., int parity = 1);
    ~Nucleus();

    void SetMass(double);
    void SetJ(double);
    void SetCharge(int);
    void SetParity(int);
    
    double Mass() const;
    double M() const;
    int Charge() const;
    int Q() const;
    double J() const;
    int Parity() const;
    int Pi() const;
    int A() const;
    int Z() const;

    /**
    * Generates a unique integer for each A,Z using the Cantor pairing function.
    * Can be used to see if two instances of the class describes the same nucleus.
    */
    int IsotopeID();
};
} // namespace rmat
#endif  //NUCLEUS_H
