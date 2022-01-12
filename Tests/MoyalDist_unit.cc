#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <chrono>
#include <getopt.h>
#include "TROOT.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TPad.h"
#include "KinKal/MatEnv/MoyalDist.hh"

void print_usage() {
  printf("Usage: MatEnv --mean f --rms f --kmax i --nsamples d \n");
}

int main(int argc, char **argv){
    double mean(3), rms(0.4);
    int kmax(10);
    long int nsamples(100000);


  static struct option long_options[] = {
    {"mean",         required_argument, 0, 'c'  },
    {"rms",          required_argument, 0, 'p'  },
    {"kmax",        required_argument, 0, 's'  },
    {"nsamples",    required_argument, 0, 'e'  },
  };

  int long_index =0;
  int opt;

  while ((opt = getopt_long_only(argc, argv,"", long_options, &long_index )) != -1) {
    switch (opt) {
      case 'c' :
        mean = atof(optarg);
        break;
      case 'p' :
        rms =atof(optarg);
        break;
      case 's' :
        kmax = atoi(optarg);
        break;
      case 'e' :
        nsamples = atol(optarg);
        break;
      default: print_usage();
               exit(EXIT_FAILURE);
    }
  }

    ////Create a root file to store the results
    std::unique_ptr<TFile> mFile( TFile::Open("MoyalDist.root", "RECREATE") );
    
    ////Create a canvas to plot the histograms and fits
    auto moyalcan = new TCanvas("moyalcan", "moyalcan", 1500, 500);
    moyalcan->Divide(3,1);

    ////Create random device for generating random number
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);

    //// Create the moyal function to fit the distribution for checking consistency
    TF1 *moyalFunc = new TF1("moyalFunc","[2] * (sqrt(1./(2.0*acos(-1.))/[1]) * exp( -0.5 * ( (x - [0])/[1] + exp(- ((x - [0])/[1]) ) )) )",0,20);
    moyalFunc->SetParameter(0,mean);
    moyalFunc->SetParameter(1,rms);
    moyalFunc->SetParameter(2,500);

    ////Create a histogram for storing distribution using Accept-Reject Method
    TH1F* histAR = new TH1F("histAR", "Accept-reject Method", 1000, 0., 20);
    MoyalDist mDist(MeanRMS(),mean, rms);
    // Start time for measuring performance
    auto start1 = std::chrono::high_resolution_clock::now(); 
    for(int i=0; i < nsamples; i++)
    {
        histAR->Fill(mDist.sampleAR());

    }
    auto finish1 = std::chrono::high_resolution_clock::now(); //end time
    std::chrono::duration<double> elapsed1 = finish1 - start1;
    std::cout << std::setw(50) << std::left << "Elapsed time for accept-reject method:" << elapsed1.count() << " s\n";
    ////Plot histAR and its fit on gPad1
    moyalcan->cd(1);
    histAR->Fit("moyalFunc","q"); 
    histAR->Draw();
    gStyle->SetOptFit(1);
    histAR->GetXaxis()->SetTitle("x");
    gPad->SetLogy();

    
    ////Create a histogram for storing distribution using InvCDF Method with default kMax=20
    TH1F* histCDFDefault = new TH1F("histCDFDefault", "InvCDF Method Default kmax(20)", 1000, 0., 20);
    // Start time for measuring performance
    auto start2 = std::chrono::high_resolution_clock::now();
    for(int i=0; i < nsamples; i++)
    {
        double sCDF = mDist.sampleInvCDF(dis(gen));
        histCDFDefault->Fill(sCDF);

    }
    auto finish2 = std::chrono::high_resolution_clock::now(); //end time
    std::chrono::duration<double> elapsed2 = finish2 - start2;
    std::cout << std::setw(50) << std::left << "Elapsed time for InvCDF default kmax(20): " << elapsed2.count() << " s\n";
    ////Plot InvCDF  and its fit on gPad1
    moyalcan->cd(2);
    histCDFDefault->Fit("moyalFunc","q"); 
    histCDFDefault->Draw();
    gStyle->SetOptFit(1);
    histCDFDefault->GetXaxis()->SetTitle("x");
    gPad->SetLogy();

    ////Create a histogram for storing distribution using InvCDF Method with user defined kmax
    TH1F* histCDFUser = new TH1F("histCDFUser", Form("InvCDF Method User kmax(%d)",kmax), 1000, 0., 20);
    MoyalDist mDistUser(MeanRMS(),mean, rms, kmax);
    // Start time for measuring performance
    auto start3 = std::chrono::high_resolution_clock::now(); 
    for(int i=0; i < nsamples; i++)
    {
         double sCDF = mDistUser.sampleInvCDF(dis(gen));
         histCDFUser->Fill(sCDF);

    }
    auto finish3 = std::chrono::high_resolution_clock::now(); //end time
    std::chrono::duration<double> elapsed3 = finish3 - start3;
    std::cout << std::setw(50) << std::left << Form("Elapsed time for InvCDF user kmax(%d): ", kmax) << elapsed3.count() << " s\n";
    ////Plot InvCDF distribution with user defined kmax and its fit on gPad1
    moyalcan->cd(3);
    histCDFUser->Fit("moyalFunc","q"); 
    histCDFUser->Draw();
    gStyle->SetOptFit(1);
    histCDFUser->GetXaxis()->SetTitle("x");
    gPad->SetLogy();
    
    histAR->Write();
    histCDFDefault->Write();
    histCDFUser->Write();
    moyalcan->Write();
 

}

