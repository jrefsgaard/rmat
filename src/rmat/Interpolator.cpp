#include <iostream>
#include "rmat/Interpolator.h"

using namespace std;

namespace rmat {

Interpolator::Interpolator(unsigned int ndata, ROOT::Math::Interpolation::Type t)
{
  type = t;
  interpolator = unique_ptr<ROOT::Math::Interpolator> (new ROOT::Math::Interpolator(ndata,t));
}

Interpolator::Interpolator(const vector<double>& vx, const vector<double>& vy, ROOT::Math::Interpolation::Type t)
{
  type = t;
  x = vx;
  y = vy;
  interpolator = SafeConstructor(vx,vy,t);
}

bool Interpolator::SetData(const vector<double>& vx, const vector<double>& vy)
{
  x = vx;
  y = vy;
  interpolator->SetData(vx,vy);
}

bool Interpolator::SetData(unsigned int ndata, const double* vx, const double* vy)
{
  x.assign(vx,vx + ndata);
  y.assign(vy,vy + ndata);
  interpolator->SetData(ndata,vx,vy);
}

Interpolator::Interpolator(const Interpolator &that)
{
  *this = that;
}

Interpolator & Interpolator::operator = (const Interpolator &that)
{
  x = that.x;
  y = that.y;
  type = that.type;
  interpolator = SafeConstructor(x,y,type);
  return *this;
}

Interpolator::Interpolator(Interpolator&& that)
{
  //To prevent code duplication we can use the move-assignment operator.
  *this = std::move(that);
}

// Move assignment operator.  
Interpolator & Interpolator::operator = (Interpolator&& that) 
{
  if(this != &that){  
    x = move(that.x);
    y = move(that.y);
    type = that.type;
    interpolator = move(that.interpolator);
    //Now 'that' is not in a valid state, so we replace its interpolator with a new one (without data).
    vector<double> x0, y0;
    that.x = x0;
    that.y = y0;
    that.interpolator = SafeConstructor(x0,y0,type);
  }
  return *this;  
} 

unique_ptr<ROOT::Math::Interpolator> Interpolator::SafeConstructor
      (const vector<double>& vx, const vector<double>& vy, ROOT::Math::Interpolation::Type t)
{
  unique_ptr<ROOT::Math::Interpolator> i;
  if(vx.size() > 0 && vy.size() > 0){
    i = unique_ptr<ROOT::Math::Interpolator> (new ROOT::Math::Interpolator(vx,vy,t));
  }
  else{
    i = unique_ptr<ROOT::Math::Interpolator> (new ROOT::Math::Interpolator(0,t));
  }
  return move(i);    
}
}//namespace rmat;
