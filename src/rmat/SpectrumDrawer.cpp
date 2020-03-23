#include <rmat/SpectrumDrawer.h>

using namespace std;

namespace rmat {

SpectrumDrawer::SpectrumDrawer(const char *name, const char *title, int nBins, int low, int high)
{
  hist = make_shared<TH1D>(name,title,nBins,low,high);
}
    
SpectrumDrawer::~SpectrumDrawer(){}
    
    
TH1D & SpectrumDrawer::Draw(Spectrum &s)
{
  int nBins = hist->GetXaxis()->GetNbins();
  double binning = hist->GetBinWidth(1);
  for(int i=1; i<=nBins; i++){
    double ei = hist->GetBinCenter(i);
    double si = s.Strength(ei);
    hist->SetBinContent(i,si * binning);
  }
  return *(hist.get());
}


TH1D & SpectrumDrawer::GetHistogram()
{
  return *(hist.get());
}
}//namespace rmat
