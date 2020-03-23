#include <vector>
#include <cmath>
#include <iostream>
#include <assert.h>
#include "rmat/Level.h"

namespace rmat {

Level::Level(double E, std::vector<double> gamma, double j, int parity)
{
  SetEnergy(E);
  widths = gamma;
  SetJ(j);
  SetParity(parity);
}

Level::~Level() {}

void Level::SetE(double E)
{
  SetEnergy(E);
}

void Level::SetEnergy(double E)
{
  energy = E;
}

void Level::SetJ(double j)
{
  jx2 = lrint(2*j);
}

void Level::SetParity(int p)
{
  assert(p == -1 || p == 1);
  pi = p;
}

double Level::E() const
{
  return energy;
}

double Level::GetEnergy() const
{
  return energy;
}

double Level::J() const
{
  return jx2/2.;
}

int Level::Parity() const
{
  return pi;
}

int Level::Pi() const
{
  return pi;
}

void Level::SetGammas(std::vector<double> gamma) 
{
  widths = gamma;
}

std::vector<double> & Level::Gammas()
{
  return widths;
}

std::vector<double> & Level::GetWidthAmplitudes()
{
  return widths;
}

void Level::SetWidthAmplitude(int c, double gammac)
{
  widths.at(c) = gammac;
}

double Level::GetWidthAmplitude(int c)
{
  return widths.at(c);
}
} //namespace rmat;
