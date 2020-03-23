#include "rmat/Nucleus.h"
#include <cmath>
#include <assert.h>
#include <ausa/constants/Mass.h>

namespace rmat {

Nucleus::Nucleus(int A, int Z, double J, int parity)
{
  a = A;
  q = Z;
  m = AUSA::Constants::isotopeMass(Z,A);
  SetJ(J);
  SetParity(parity);
}

Nucleus::~Nucleus(){}

void Nucleus::SetMass(double mass)
{
  m = mass;
}

void Nucleus::SetCharge(int Q)
{
  q = Q;
}

void Nucleus::SetJ(double j)
{
  jx2 = lrint(2.* j);
}

void Nucleus::SetParity(int parity)
{
  assert(parity == -1 || parity == 1);
  pi = parity;
}

double Nucleus::Mass() const
{
  return m;
}

double Nucleus::M() const
{
  return m;
}

int Nucleus::Charge() const
{
  return q;
}

int Nucleus::Q() const
{
  return q;
}

double Nucleus::J() const
{
  return jx2/2.;;
}

int Nucleus::Parity() const
{
  return pi;
}

int Nucleus::Pi() const
{
  return pi;
}

int Nucleus::A() const
{
  return a;
}

int Nucleus::Z() const
{
  return q;
}

int Nucleus::IsotopeID()
{
  return lrint(0.5 * (a + q) * (a + q + 1) + q);
}
} //namespace rmat
