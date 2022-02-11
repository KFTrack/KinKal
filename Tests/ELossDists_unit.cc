//
// test basic functions of MatEnv and materials
//
#include "KinKal/MatEnv/MatDBInfo.hh"
#include "KinKal/MatEnv/DetMaterial.hh"
#include "KinKal/MatEnv/SimpleFileFinder.hh"

#include <iostream>
#include <stdio.h>
#include <iostream>
#include <getopt.h>
#include <iomanip>
#include <random>
#include <cmath>
#include <chrono>

#include "TH1F.h"
#include "TSystem.h"
#include "TFile.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "KinKal/MatEnv/ELossDistributions.hh"

using namespace std;
using namespace MatEnv;
using KinKal::MoyalDist;
using KinKal::BremssLoss;

void print_usage() {
  printf("Usage: MatEnv --material c --particle i --momentum f  --thickness f\n");
}

int main(int argc, char **argv) {

  string matname("straw-wall");
  double momentum(100.0);
  double thickness(0.0015);
  int imass(0);
  double masses[5]={0.511,105.66,139.57, 493.68, 938.0};
  const char* pnames[5] = {"electron","muon","pion","kaon","proton"};
  double pmass;
  string pname;

  static struct option long_options[] = {
    {"material",     required_argument, 0, 'c'  },
    {"particle",     required_argument, 0, 'p'  },
    {"momentum",     required_argument, 0, 's'  },
    {"thickness",     required_argument, 0, 't'  },
  };

  int long_index =0;
  int opt;

  while ((opt = getopt_long_only(argc, argv,"", long_options, &long_index )) != -1) {
    switch (opt) {
      case 'c' :
        matname = string(optarg);
        break;
      case 'p' :
        imass =atoi(optarg);
        break;
      case 's' :
        momentum = atof(optarg);
        break;
      case 't' :
        thickness = atof(optarg);
        break;
      default: print_usage();
               exit(EXIT_FAILURE);
    }
  }
  pmass = masses[imass];
  pname = pnames[imass];
  cout << "Test for particle " << pname  << " mass " << pmass << endl;
  cout << "Searching for material " << matname << endl;
  MatEnv::SimpleFileFinder sfinder;
  MatDBInfo matdbinfo(sfinder,MatEnv::DetMaterial::moyalmean);
  const DetMaterial* dmat = matdbinfo.findDetMaterial(matname);
  std::cout << "Found DetMaterial " << dmat->name() << std::endl;
  
  double eloss = dmat->energyLoss(momentum,thickness,pmass);
  double elossrms = dmat->energyLossRMS(momentum,thickness,pmass);
  MoyalDist mDist(MoyalDist::MeanRMS(abs(eloss), elossrms));
  std::cout << "abs(eloss) = " << abs(eloss)  << std::endl;
  std::cout << "elossrms = " << elossrms  << std::endl;
  BremssLoss bLoss;

  double radFrac = dmat->radiationFraction(thickness);
  std::cout << "radiation fraction == " << radFrac << std::endl;

  
  std::unique_ptr<TFile> mFile( TFile::Open("ELossDists.root", "RECREATE") );
  TH1F* histBrem = new TH1F("histBrem", "Bremss Loss", 500, 0, 100e-3);
  TH1F* histCol = new TH1F("histCol", "Collision Loss", 500,  0., 100e-3);
  TH1F* elossTotal = new TH1F("elossTotal", "Total Loss", 500,  0., 100e-3);
  
  double nSamples = 100000;

  for (int i = 0; i < nSamples; i++){

    // Here we assume that the particle loses energy first in collision and then in radiation
    // Changing the order doesn't change the results. 
    double colloss = mDist.sampleAR();
    double bremloss = bLoss.sampleSSPGamma(momentum - colloss,radFrac); 

    histBrem->Fill(bremloss);
    histCol->Fill(colloss);
    elossTotal->Fill(bremloss+colloss);
  } 

  histBrem->SetLineColor(kRed);
  histBrem->SetFillColor(kRed);
  histBrem->SetFillStyle(3001);
  histBrem->Write();

  histCol->SetLineColor(kBlue);
  histCol->SetFillColor(kBlue);
  histCol->SetFillStyle(3001);
  histCol->Write();

  elossTotal->SetLineColor(kGreen);
  elossTotal->SetFillColor(kGreen);
  elossTotal->SetFillStyle(3001);
  elossTotal->Write();

  TCanvas* canvas = new TCanvas("canvas","canvas",1000,1000);
  canvas->cd();
  elossTotal->Draw();
  
  histCol->Draw("same");
  histBrem->Draw("same");
  gPad->SetLogy();
  gPad->BuildLegend();
  canvas->Write();

  // Uncheck the lines below to compare performance of two BremssLoss methods
  // //Check which method generates faster 
  // auto start1 = std::chrono::high_resolution_clock::now(); //start time
  // for (int i = 0; i < nSamples; i++){
  //   bLoss.sampleSTDGamma(momentum,radFrac);
  // }
  // auto finish1 = std::chrono::high_resolution_clock::now(); //end time
  // std::chrono::duration<double> elapsed1 = finish1 - start1;
  // std::cout << std::setw(50) << std::left << "Elapsed time for STDGamma method:" << elapsed1.count() << " s\n";

  // auto start2 = std::chrono::high_resolution_clock::now(); //start time
  // for (int i = 0; i < nSamples; i++){
  //   bLoss.sampleSSPGamma(momentum,radFrac);
  // }
  // auto finish2 = std::chrono::high_resolution_clock::now(); //end time
  // std::chrono::duration<double> elapsed2 = finish2 - start2;
  // std::cout << std::setw(50) << std::left << "Elapsed time for SSPGamma method:" << elapsed2.count() << " s\n";

}
