#include <rmat/threebody/ObservedDalitzPlot.h>

using namespace std;

namespace rmat::threebody {

ObservedDalitzPlot::ObservedDalitzPlot(shared_ptr<DalitzPlot> model_,
                                       shared_ptr<TH3D> simulated,
                                       shared_ptr<TH3D> accepted)
: model(model_)
{
  if((simulated->GetNbinsX() != accepted->GetNbinsX()) ||
     (simulated->GetNbinsY() != accepted->GetNbinsY()) ||
     (simulated->GetNbinsZ() != accepted->GetNbinsZ())){
    cout << "ObservedDalitzPlot::ObservedDalitzPlot():" ;
    cout << "Binning of response histograms doesn't match." << endl;
    exit(EXIT_FAILURE);
  }

  acceptance = shared_ptr<TH3D>((TH3D*)accepted->Clone());
  acceptance->Divide(simulated.get());
}

void ObservedDalitzPlot::SetModel(shared_ptr<DalitzPlot> model_)
{
  model = model_;
}

DalitzPlot & ObservedDalitzPlot::GetModel()
{
  return *(model.get());
}

shared_ptr<TH3D> ObservedDalitzPlot::GetSpectrum()
{
  shared_ptr<TH3D> result((TH3D*)acceptance->Clone());
  for(int binx=1; binx<=result->GetXaxis()->GetNbins(); binx++){
    for(int biny=1; biny<=result->GetYaxis()->GetNbins(); biny++){
      for(int binz=1; binz<=result->GetZaxis()->GetNbins(); binz++){
        double acc = result->GetBinContent(binx,biny,binz);
        double x = result->GetXaxis()->GetBinCenter(binx);
        double y = result->GetYaxis()->GetBinCenter(biny);
        double Q = result->GetZaxis()->GetBinCenter(binz);
        double f = model->Value(x,y,Q);
        result->SetBinContent(binx,biny,binz,f*acc);
      }
    }
  }
  return result;
}

void ObservedDalitzPlot::SetParameters(std::vector<double> par)
{
  model->SetParameters(par);
}

void ObservedDalitzPlot::PrintParameters()
{
  model->GetDecayWeight().PrintParameters();
}

int ObservedDalitzPlot::NDim()
{
  return model->NDim();
}
}
