#include <TMath.h>
#include <rmat/LikelihoodChiSquare.h>
#include <rmat/RMatrixUtils.h>

using namespace std;
using namespace TMath;

namespace rmat {

LikelihoodChiSquare::LikelihoodChiSquare(shared_ptr<TH1> d, std::shared_ptr<Spectrum> m) : ChiSquare(move(d),move(m)) {}

LikelihoodChiSquare::~LikelihoodChiSquare(){}

double LikelihoodChiSquare::Evaluate() const
{
  nBins = 0;
  double sum = 0.;
  double counts = data->Integral(minBin,maxBin);
  for(int i=minBin; i<=maxBin; i++){
    double xi = data->GetBinCenter(i);
    double yi = binning * model->Strength(xi);
    double ni = data->GetBinContent(i);
    if(IsAlmostZero(yi)) continue;
    nBins ++;
    if(ni <= 0) sum += yi;
    else sum += ni * Log(ni / yi) + yi - ni;
  }
  return 2 * sum;
}

ROOT::Math::IBaseFunctionMultiDim * LikelihoodChiSquare::Clone() const
{
  LikelihoodChiSquare *clone = new LikelihoodChiSquare(*this);  
  return clone;
}
}//namespace rmat;
