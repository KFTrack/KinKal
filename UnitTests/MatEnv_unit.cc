// 
// test basic functions of MatEnv and materials
//
#include "MatEnv/MatDBInfo.hh"
#include "MatEnv/DetMaterial.hh"

#include <iostream>
#include <stdio.h>
#include <iostream>
#include <getopt.h>

#include "TH1F.h"
#include "TSystem.h"
#include "THelix.h"
#include "TFile.h"
#include "TPolyLine3D.h"
#include "TAxis3D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TPolyLine3D.h"
#include "TPolyMarker3D.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TProfile2D.h"

using namespace std;
using namespace MatEnv;

void print_usage() {
  printf("Usage: MatEnv --material c --particle i --momstart f --momend f --thickness f\n");
}

int main(int argc, char **argv) {

  string matname("straw-wall");
  double momstart(10.0), momend(200.0);
  double thickness(0.015);
  int imass(0);
  double masses[5]={0.511,105.66,139.57, 493.68, 938.0};
  const char* pnames[5] = {"electron","muon","pion","kaon","proton"};
  double pmass;
  string pname;

  static struct option long_options[] = {
    {"material",     required_argument, 0, 'c'  },
    {"particle",     required_argument, 0, 'p'  },
    {"momstart",     required_argument, 0, 's'  },
    {"momend",     required_argument, 0, 'e'  },
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
	momstart = atof(optarg);
	break;
      case 'e' : 
	momend = atof(optarg);
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
  MatDBInfo matdbinfo;
  const DetMaterial* dmat = matdbinfo.findDetMaterial(matname);
  if(dmat != 0){
    cout << "Found DetMaterial " << dmat->name() << endl;
    unsigned nstep(100);
    double momstep = (momend-momstart)/(nstep-1); 
    TGraph* geloss = new TGraph(nstep);
    string title = string("Eloss vs Momentum ")
     + dmat->name() + string(" ") + pname
     + string(";Mom (MeV/c);MeV");
    geloss->SetTitle(title.c_str());
    TGraph* gelossrms = new TGraph(nstep);
    title = string("Eloss RMS vs Momentum ")
     + dmat->name() + string(" ") + pname
     + string(";Mom (MeV/c);MeV");
    gelossrms->SetTitle(title.c_str());
    TGraph* gascat = new TGraph(nstep);
    title = string("Scattering RMS vs Momentum ") 
     + dmat->name() + string(" ") + pname
     + string(";Mom (MeV/c);Radians");
    gascat->SetTitle(title.c_str());
    TGraph* gbetagamma = new TGraph(nstep);
    title = string("Particle #beta#gamma vs Momentum ") 
     + dmat->name() + string(" ") + pname
     + string(";Mom (MeV/c);#beta#gamma");
    gbetagamma->SetTitle(title.c_str());
    for(unsigned istep = 0;istep < nstep; istep++){
      double mom = momstart + istep*momstep;
      double eloss = dmat->energyLoss(mom,thickness,pmass);
      geloss->SetPoint(istep,mom,eloss);
      double elossrms = dmat->energyLossRMS(mom,thickness,pmass);
      gelossrms->SetPoint(istep,mom,elossrms);
      double ascat = dmat->scatterAngleRMS(mom,thickness,pmass);
      gascat->SetPoint(istep,mom,ascat);
      double betagamma = dmat->particleBetaGamma(mom,pmass);
      gbetagamma->SetPoint(istep,mom,betagamma);
    }
    TFile mefile("MatEnv.root","RECREATE");
    TCanvas* matcan = new TCanvas("matcan","MatEnv",1000,1000);
    matcan->Divide(2,2);
    matcan->cd(1);
    geloss->Draw("AC*");
    matcan->cd(2);
    gelossrms->Draw("AC*");
    matcan->cd(3);
    gascat->Draw("AC*");
    matcan->cd(4);
    gbetagamma->Draw("AC*");
    matcan->Write();
    mefile.Write();
    mefile.Close();
  }
}
