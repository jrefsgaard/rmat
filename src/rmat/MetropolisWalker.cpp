#include "rmat/MetropolisWalker.h"
#include <iostream>
#include <limits>
#include <stdlib.h>
#include <cmath>
#include <TRandom3.h>
#include <TTree.h>
#include <TFile.h>

using namespace std;

namespace rmat {

MetropolisWalker::MetropolisWalker(ROOT::Math::Functor f, int scaling)
: chi2(f), Nscale(scaling)
{
  T = 1.;
  int dim = chi2.NDim();
  ParameterInfo defaultOption;
  defaultOption.fixed = false;
  defaultOption.sigma = 1.;
  defaultOption.min = -numeric_limits<double>::max();
  defaultOption.max = numeric_limits<double>::max();
  options = vector<ParameterInfo>(dim,defaultOption);
}

MetropolisWalker::~MetropolisWalker() {}

double MetropolisWalker::Evaluate(vector<double> x)
{
  double val = chi2(x.data());
  Step step;
  step.x = x;
  step.value = val;
  history.push_back(step);
  cout << val ;
  for(double xi : x){ cout << "  " << xi ;}
  cout << endl;
  return val;
}

void MetropolisWalker::Walk(vector<double> x0, int steps)
{
  int dim = chi2.NDim();
  //Check for valid starting point.
  if(x0.size() != dim){
    cout << "MetropolisWalker::Walk(): x0 has wrong dimension, size = " << x0.size() << "!" << endl;
    exit(EXIT_FAILURE);
  }
  
  //Evaluate and store the starting point.
  double f0 = Evaluate(x0);
  Step step0;
  step0.x = x0;
  step0.value = f0;
  position.push_back(step0);
  
  //We keep track of the acceptance ratio for each parameter.
  vector<int> accepted (dim,0);
  
  //Start the random walk.
  for(int i=0; i<steps; i++){
    position.push_back(position.back()); //Copy the (i-1)'th entry.
    vector<double> & xi = position.back().x;
    for(int j=0; j<dim; j++){
      ParameterInfo & option = options.at(j);
      if(option.fixed) continue;  //If fixed, go on to next parameter.
      TRandom3 rGen(0);
      double step = rGen.Gaus(0.,option.sigma);
      xi.at(j) += step;
      if(xi.at(j) < option.min || xi.at(j) > option.max){
        //If outside boundary, reject the step an go on to next parameter.
        xi.at(j) -= step;
        continue;
      }
      //Ok, the proposed step is valid. We evaluate the function...
      double fOld = position.back().value;
      double fNew = Evaluate(xi);
      double p = fmin(1.,exp(-(fNew - fOld)/T)); //Probability of taking the step
            
      //Make a decision.
      if(rGen.Uniform() < p){
        //Complete the step.
        position.back().value = fNew; 
        accepted.at(j) += 1;
      }
      else{
         //Revert the step.
        xi.at(j) -= step;
      }
    }
    
    //Every Nscale steps we check the acceptance rate and modify the step size accordingly.
    if(i % Nscale == 0){
      for(int j=0; j<dim; j++){
        int Nacc = accepted.at(j);
        ParameterInfo & opt = options.at(j);
        //cout << opt.sigma ;
        if(opt.fixed) continue;
        if(Nacc == 0){
          opt.sigma = opt.sigma/3.;  //A bit arbitrary.
        }
        else{
          double accRate = (double)Nacc/Nscale;
          //cout << "  " << accRate ;
          if(accRate < 0.15 || accRate > 0.40){
            opt.sigma = opt.sigma * accRate / 0.25;
          }
        }
        //cout << endl;
        accepted.at(j) = 0;
      }
    }      
  }
}
    
void MetropolisWalker::SetFunction(ROOT::Math::Functor f)
{
  chi2 = f;
}

void MetropolisWalker::SetTemperature(double temp)
{
  T = temp;
}

void MetropolisWalker::SetScalingPeriod(int nsteps)
{
  Nscale = nsteps;
}

void MetropolisWalker::SetParOptions(int ipar, double sigma, double pmin, double pmax)
{
  ParameterInfo & opt = options.at(ipar);
  opt.sigma = sigma;
  opt.min = pmin;
  opt.max = pmax;
}

void MetropolisWalker::SetParOptions(vector<ParameterInfo> opts)
{
  options = opts;
}

void MetropolisWalker::FixParameter(int ipar)
{
  options.at(ipar).fixed = true;
}

void MetropolisWalker::ReleasePar(int ipar)
{
  options.at(ipar).fixed = false;
}
    
ROOT::Math::Functor & MetropolisWalker::GetFunction()
{
  return chi2;
}

double MetropolisWalker::GetTemperature()
{
  return T;
}

int MetropolisWalker::GetScalingPeriod()
{
  return Nscale;
}

vector<Step> & MetropolisWalker::GetSteps()
{
  return position;
}

vector<Step> & MetropolisWalker::GetHistory()
{
  return history;
}

ParameterInfo & MetropolisWalker::GetParOptions(int ipar)
{
  return options.at(ipar);
}

vector<ParameterInfo> & MetropolisWalker::GetParOptions()
{
  return options;
}

void MetropolisWalker::SaveTree(string fileName)
{
  TFile outFile(fileName.c_str(),"RECREATE");
  TTree *tree = new TTree("walk","Result of Metropolis algorithm.");
  int dim = position.at(0).x.size();
  double val;
  tree->Branch("val",&val);
  double *p = new double(dim);
  //cout << "Creating branches..." << endl;
  for(int i=0; i<dim; i++){
    char branchName[32];
    sprintf(branchName,"p%d",i);
    tree->Branch(branchName,p+i);
  }
  //cout << "Filling the tree..." << endl;
  for(Step & step : position){
    val = step.value;
    for(int i=0; i<dim; i++){
      p[i] = step.x[i];
    }
    tree->Fill();
  }
  //cout << "Writing to file..." << endl;
  tree->Write();
  outFile.Close();
}
}
