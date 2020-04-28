// 
// test basic functions of BField class
//
#include "KinKal/LHelix.hh"
#include "KinKal/BField.hh"

#include <iostream>
#include <stdio.h>
#include <iostream>
#include <getopt.h>

#include "TFile.h"
#include "TH1F.h"
#include "TSystem.h"
#include "THelix.h"
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

using namespace KinKal;
using namespace std;

void print_usage() {
  printf("Usage: BField  --momentum f --costheta f --phi f --charge i --bz f --dbz f --dby f --bgrad f --tmin f --tmax f \n");
}

int main(int argc, char **argv) {
  typedef LHelix KTRAJ;
  double mom(105.0), cost(0.7), phi(0.5);
  int icharge(-1);
  double pmass(0.511);
  double tmin(-10.0), tmax(10.0);
  BField *BF(0);
  float Bgrad(0.0), dBx(0.0), dBy(0.0), dBz(0.0);
  static struct option long_options[] = {
    {"momentum",     required_argument, 0, 'm' },
    {"costheta",     required_argument, 0, 'c'  },
    {"phi",     required_argument, 0, 'f'  },
    {"charge",     required_argument, 0, 'q'  },
    {"dBx",     required_argument, 0, 'x'  },
    {"dBy",     required_argument, 0, 'y'  },
    {"dBz",     required_argument, 0, 'Z'  },
    {"Bgrad",     required_argument, 0, 'g'  },
    {"tmin",     required_argument, 0, 't'  },
    {"tmax",     required_argument, 0, 'd'  },
  };
  int opt;
  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"",
	  long_options, &long_index )) != -1) {
    switch (opt) {
      case 'm' : mom = atof(optarg);
		 break;
      case 'c' : cost = atof(optarg);
		 break;
      case 'f' : phi = atof(optarg);
		 break;
      case 'x' : dBx = atof(optarg);
		 break;
      case 'y' : dBy = atof(optarg);
		 break;
      case 'Z' : dBz = atof(optarg);
		 break;
      case 'g' : Bgrad = atof(optarg);
		 break;
      case 't' : tmin = atof(optarg);
		 break;
      case 'd' : tmax = atof(optarg);
		 break;
      default: print_usage();
	       exit(EXIT_FAILURE);
    }
  }
  Vec3 bnom(0.0,0.0,1.0);
  Vec3 bsim;
  if(Bgrad != 0){
  } else {
    Vec3 bsim(dBx,dBy,1.0+dBz);

    BF = new UniformBField(bsim);
    bnom = Vec3(0.0,0.0,1.0);
  }
  Vec4 origin(0.0,0.0,0.0,0.0);
  float sint = sqrt(1.0-cost*cost);
  Mom4 momv(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost,pmass);
  // first, create a traj based on the actual field at this point
  KTRAJ simt(origin,momv,icharge,bsim);
  // then, a traj based on the nominal field at this point
  KTRAJ nomt(origin,momv,icharge,bnom);


  TFile tpfile("BField.root","RECREATE");

  return 0;
}


