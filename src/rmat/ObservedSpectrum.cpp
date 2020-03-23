#include <rmat/ObservedSpectrum.h>
#include <iostream>

using namespace std;
using namespace arma;

namespace rmat {

ObservedSpectrum::ObservedSpectrum(shared_ptr<Spectrum> m, TMatrixD &r, double Emin, double Emax)
: model(m), emin(Emin), emax(Emax)
{ 
  rebin = 1;
  SetResponse(r);
}

ObservedSpectrum::~ObservedSpectrum() {}

//The working horse of this class.
void ObservedSpectrum::CalculateSpectrum()
{
  int nbins = response.n_cols;

  //First, we update the theoretical spectrum.
  raw.set_size(nbins);
  for(int i=0; i<nbins; i++){
    double ei = Bin2Energy(i);
    double si = model->Strength(ei);
    //cout << "Calc: i = " << i << ",  ei = " << ei << ",  si = " << si << endl;
    raw(i) = si;
  }
  
  //Next, we apply the response matrix.
  
  //cout << "nrows = " << response.n_rows << ",  ncols = " << response.n_cols ;
  //cout << ",  raw.n_rows = " << raw.n_rows << endl;
  
  observed = response * raw;
  
  //cout << ",  observed.n_rows = " << observed.n_rows << endl;
  
  RebinSpectrum();
}

//Rebins the observed spectrum.
void ObservedSpectrum::RebinSpectrum()
{
  if(rebin == 1) return;
  int finalSize = (int)floor(observed.n_rows / rebin);
  //cout << " n_cols = " << observed.n_rows << ",  rebin = " << rebin << ",  finalSize = " << finalSize << endl;
  Col<double> result(finalSize);
  for(int i=0; i<finalSize; i++){
    int start = i * rebin;
    int end = (i+1) * rebin - 1;
    //cout << "start = " << start << ",  end = " << end << endl;
    double content = sum(observed.subvec(start,end));
    double mean = content / rebin;   //Preserve the 'height' of the spectrum.
    result(i) = mean;
  }
  observed = result;
}

void ObservedSpectrum::SetResponse(TMatrixD &r)
{
  int nrows = r.GetNrows();
  int ncols = r.GetNcols();
  //cout << "SetResponse() rows = " << nrows << ",  cols = " << ncols << endl;
  if(nrows == 0 || ncols == 0){
    cout << "  ObservedSpectrum::SetResponse(): Dimension issue, nrows = " << nrows << ", ncols = " << ncols << endl;
    exit(EXIT_FAILURE);
  }   
  response.set_size(nrows,ncols);
  for(int i=0; i<nrows; i++){
    for(int j=0; j<ncols; j++){
      response(i,j) = r(i,j);  //element-wise copy. Probably very slow...?
    }
  }
  
  //Response has changed, update spectrum.
  CalculateSpectrum();
}

//Returns the center of the energy bin 'i'.
double ObservedSpectrum::Bin2Energy(int i)
{
  int nbins = response.n_cols;
  if(i < 0 || i >= nbins){
    cout << "  ObservedSpectrum::Bin2Energy(): Bin does not exist, i = " << i << ", nbins = " << nbins << endl;
    exit(EXIT_FAILURE);
  }  
  return (emax - emin) / nbins * (i + 0.5) + emin;
}

int ObservedSpectrum::Energy2Bin(double E)
{
  if(E < emin || E > emax){
    cout << "  ObservedSpectrum::Energy2Bin(): Energy " << E << " out of range!" << endl;
    exit(EXIT_FAILURE);
  }
  
  int nbins = observed.n_rows;
  //cout << "E2B: nbins = " << nbins << ",  emin = " << emin << ",  emax = " << emax << endl;
  return lrint(nbins*(E-emin)/(emax - emin)-0.5);
}

double ObservedSpectrum::Strength(double Ec)
{
  int i = Energy2Bin(Ec);
  //cout << " Strength: Ec = " << Ec << ",  i = " << i << endl;
  double si = observed(i);
  return si;
}

void ObservedSpectrum::SetParameters(std::vector<double> par)
{
  model->SetParameters(par);
  CalculateSpectrum();  //The parameters have changed, so we must update all spectra.
}

void ObservedSpectrum::PrintParameters()
{
  model->PrintParameters();
}

int ObservedSpectrum::NDim()
{
  return model->NDim();
}
    
void ObservedSpectrum::Rebin(int N)
{
  rebin = N;
  RebinSpectrum();
}

}//namespace rmat
