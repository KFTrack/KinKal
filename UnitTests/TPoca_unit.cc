// 
// test basic functions of TPoca using LHelix and TLine
//
#include "KinKal/LHelix.hh"
#include "KinKal/TLine.hh"
#include "KinKal/TPoca.hh"
#include "KinKal/BField.hh"
#include "CLHEP/Units/PhysicalConstants.h"

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

using namespace KinKal;
using namespace std;
// avoid confusion with root
using KinKal::TLine;

void print_usage() {
  printf("Usage: TPoca --charge i--gap f --time f --dtime f --dphi f --vprop f\n");
}

int main(int argc, char **argv) {
  int opt;
  double mom(105.0), cost(0.7), phi(0.5);
  int icharge(-1);
  double pmass(0.511), oz(100.0), ot(0.0);
  double time(8.0), dtime(20.0), dphi(5.0);
  double hlen(500.0); // half-length of the wire
  double gap(2.0); // distance between TLine and LHelix 
  double vprop(0.7);
  unsigned ndtest(50);

  static struct option long_options[] = {
    {"charge",     required_argument, 0, 'q'  },
    {"gap",     required_argument, 0, 'g'  },
    {"time",     required_argument, 0, 't'  },
    {"dtime",     required_argument, 0, 'd'  },
    {"dphi",     required_argument, 0, 'f'  },
    {"vprop",     required_argument, 0, 'v'  }
  };

  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"", 
	  long_options, &long_index )) != -1) {
    switch (opt) {
      case 'q' : icharge = atoi(optarg);
		 break;
      case 'g' : gap = atof(optarg); 
		 break;
      case 't' : time = atof(optarg);
		 break;
      case 'd' : dtime = atof(optarg);
		 break;
      case 'f' : dphi = atof(optarg);
		 break;
      case 'v' : vprop = atof(optarg);
		 break;
      default: print_usage(); 
	       exit(EXIT_FAILURE);
    }
  }
// create helix
  Vec3 bnom(0.0,0.0,1.0);
  UniformBField BF(bnom); // 1 Tesla
  Vec4 origin(0.0,0.0,oz,ot);
  float sint = sqrt(1.0-cost*cost);
  Mom4 momv(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost,pmass);
  LHelix lhel(origin,momv,icharge,bnom);
// create tline perp to z axis at the specified time, separated by the specified gap
  Vec3 pos, dir;
  lhel.position(time,pos);
  lhel.direction(time,dir);
  // rotate the direction
  double lhphi = atan2(dir.Y(),dir.X());
  double pphi = lhphi + dphi;
  Vec3 pdir(cos(pphi),sin(pphi),0.0);
  double pspeed = CLHEP::c_light*vprop; // vprop is relative to c
  Vec3 pvel = pdir*pspeed;
  // shift the position
  Vec3 perpdir(-sin(phi),cos(phi),0.0);
  Vec3 ppos = pos + gap*perpdir;
// time range;
  double ptime = time+dtime;
  TRange prange(ptime-hlen/pspeed, ptime+hlen/pspeed);
  // create the TLine
  TLine tline(ppos, pvel,time+dtime,prange);
  // create TPoca from these
  TPoca<LHelix,TLine> tp(lhel,tline);
  cout << "TPoca status " << tp.statusName() << " doca " << tp.doca() << " dt " << tp.deltaT() << endl;
  Vec3 thpos, tlpos;
  tp.particleTraj().position(tp.particlePoca().T(),thpos);
  tp.sensorTraj().position(tp.sensorPoca().T(),tlpos);
  double refd = tp.doca();
  cout << " Helix Pos " << pos << " TPoca LHelix pos " << thpos << " TPoca TLine pos " << tlpos << endl;
  cout << " TPoca particlePoca " << tp.particlePoca() << " TPoca sensorPoca " << tp.sensorPoca()  << " DOCA " << refd << endl;

  // now derivatives
  TPoca<LHelix,TLine> tdp(tp);
  cout << "TPoca dDdP " << tdp.dDdP() << " dTdP " << tdp.dTdP() << endl;
  // test against numerical derivatives
  std::vector<TGraph*> dtpoca;
  // range to change specific parameters; most are a few mm
  std::vector<double> pchange = {10.0,5.0,10.0,10.0,0.1,0.1};
  TFile tpfile("TPoca.root","RECREATE");
  TCanvas* dtpcan = new TCanvas("dtpcan","DTPoca",1200,800);
  dtpcan->Divide(3,2);
  for(int ipar=0;ipar<lhel.npars_;ipar++){
    LHelix::ParamIndex parindex = static_cast<LHelix::ParamIndex>(ipar);
    dtpoca.push_back(new TGraph(ndtest));
    string ts = LHelix::paramTitle(parindex)+string(" DOCA Change;#Delta DOCA (exact);#Delta DOCA (derivative)");
    dtpoca.back()->SetTitle(ts.c_str());
    double dstep = pchange[ipar]/(ndtest-1);
    double dstart = -0.5*pchange[ipar];
    for(unsigned istep=0;istep<ndtest;istep++){
   // compute exact change in DOCA 
      auto dvec = lhel.params().parameters();
      double dpar = dstart + dstep*istep;
      dvec[ipar] += dpar; 
      LHelix::PDATA pdata(dvec,lhel.params().covariance());
      LHelix dlhel(pdata,lhel.mass(),lhel.charge(),bnom);
      TPoca<LHelix,TLine> dtp(dlhel,tline);
      double xd = dtp.doca();
      // now derivatives
      double dd = tdp.dDdP()[ipar]*dpar;
      dtpoca.back()->SetPoint(istep,xd-refd,dd);
    }
    dtpcan->cd(ipar+1);
    dtpoca.back()->Draw("AC*");
  }
  dtpcan->Write();
  tpfile.Write();
  tpfile.Close();
  return 0;
}


