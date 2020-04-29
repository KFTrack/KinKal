// 
// test basic functions of TPoca using KTraj and TLine
//
#include "KinKal/LHelix.hh"
#include "KinKal/IPHelix.hh"
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
  printf("Usage: TPoca --charge i--gap f --tmin f --tmax f --phi f --vprop f --costheta f --eta f\n");
}

int main(int argc, char **argv) {
  //  typedef IPHelix KTRAJ;
  typedef LHelix KTRAJ;
  typedef TPoca<KTRAJ,TLine> TPOCA;
  int opt;
  double mom(105.0), cost(0.7), phi(0.5);
  int icharge(-1);
  double pmass(0.511), oz(0.0), ot(0.0);
  double tmin(-10.0), tmax(10.0);
  double hlen(500.0); // half-length of the wire
  double gap(2.0); // distance between TLine and KTRAJ 
  double vprop(0.7);
  double eta(0.0);
  unsigned nstep(50),ntstep(10);

  static struct option long_options[] = {
    {"charge",     required_argument, 0, 'q'  },
    {"costheta",     required_argument, 0, 'c'  },
    {"gap",     required_argument, 0, 'g'  },
    {"tmin",     required_argument, 0, 't'  },
    {"tmax",     required_argument, 0, 'd'  },
    {"phi",     required_argument, 0, 'f'  },
    {"vprop",     required_argument, 0, 'v'  },
    {"eta",     required_argument, 0, 'e'  },
  };

  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"", 
	  long_options, &long_index )) != -1) {
    switch (opt) {
      case 'c' : cost = atof(optarg);
		 break;
      case 'q' : icharge = atoi(optarg);
		 break;
      case 'g' : gap = atof(optarg); 
		 break;
      case 't' : tmin = atof(optarg);
		 break;
      case 'd' : tmax = atof(optarg);
		 break;
      case 'f' : phi = atof(optarg);
		 break;
      case 'v' : vprop = atof(optarg);
		 break;
      case 'e' : eta = atof(optarg);
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
  KTRAJ lhel(origin,momv,icharge,bnom);
  TFile tpfile("TPoca.root","RECREATE");
  TCanvas* ttpcan = new TCanvas("ttpcan","DToca",1200,800);
  ttpcan->Divide(3,2);
  TCanvas* dtpcan = new TCanvas("dtpcan","DDoca",1200,800);
  dtpcan->Divide(3,2);
  std::vector<TGraph*> dtpoca, ttpoca;
  std::vector<double> pchange = {10.0,5.0,10.0,10.0,0.1,0.1};
  for(int ipar=0;ipar<lhel.npars_;ipar++){
    KTRAJ::ParamIndex parindex = static_cast<KTRAJ::ParamIndex>(ipar);
    dtpoca.push_back(new TGraph(nstep*ntstep));
    string ts = KTRAJ::paramTitle(parindex)+string(" DOCA Change;#Delta DOCA (exact);#Delta DOCA (derivative)");
    dtpoca.back()->SetTitle(ts.c_str());
    ttpoca.push_back(new TGraph(nstep*ntstep));
    ts = KTRAJ::paramTitle(parindex)+string(" TOCA Change;#Delta TOCA (exact);#Delta TOCA (derivative)");
    ttpoca.back()->SetTitle(ts.c_str());
  }
  for(unsigned itime=0;itime < ntstep;itime++){
    float time = tmin + itime*(tmax-tmin)/(ntstep-1);
    // create tline perp to trajectory at the specified time, separated by the specified gap
    Vec3 pos, dir, perp1, perp2;
    lhel.position(time,pos);
    lhel.direction(time,dir);
    KTRAJ::DVEC der;
    lhel.momDeriv(KInter::theta1,time,der,perp1); 
    lhel.momDeriv(KInter::theta2,time,der,perp2); 
    // choose a specific direction for DOCA
    // the line traj must be perp. to this and perp to the track
    Vec3 docadir = cos(eta)*perp1 + sin(eta)*perp2;
    Vec3 pdir = sin(eta)*perp1 - cos(eta)*perp2;
    double pspeed = CLHEP::c_light*vprop; // vprop is relative to c
    Vec3 pvel = pdir*pspeed;
    // shift the position
    Vec3 ppos = pos + gap*docadir;
    // time range;
    TRange prange(time-hlen/pspeed, time+hlen/pspeed);
    // create the TLine
    TLine tline(ppos, pvel,time,prange);
    // create TPoca from these
    TPOCA tp(lhel,tline);
//    cout << "TPoca status " << tp.statusName() << " doca " << tp.doca() << " dt " << tp.deltaT() << endl;
    Vec3 thpos, tlpos;
    tp.particleTraj().position(tp.particlePoca().T(),thpos);
    tp.sensorTraj().position(tp.sensorPoca().T(),tlpos);
    double refd = tp.doca();
    double reft = tp.deltaT();
//    cout << " Helix Pos " << pos << " TPoca KTRAJ pos " << thpos << " TPoca TLine pos " << tlpos << endl;
//    cout << " TPoca particlePoca " << tp.particlePoca() << " TPoca sensorPoca " << tp.sensorPoca()  << " DOCA " << refd << endl;
//    cout << "TPoca dDdP " << tp.dDdP() << " dTdP " << tp.dTdP() << endl;
    // test against numerical derivatives
    // range to change specific parameters; most are a few mm
    for(int ipar=0;ipar<lhel.npars_;ipar++){
      double dstep = pchange[ipar]/(nstep-1);
      double dstart = -0.5*pchange[ipar];
      for(unsigned istep=0;istep<nstep;istep++){
	// compute exact change in DOCA 
	auto dvec = lhel.params().parameters();
	double dpar = dstart + dstep*istep;
	dvec[ipar] += dpar; 
	KTRAJ::PDATA pdata(dvec,lhel.params().covariance());
	KTRAJ dlhel(pdata,lhel.mass(),lhel.charge(),bnom);
	TPOCA dtp(dlhel,tline);
	double xd = dtp.doca();
	double xt = dtp.deltaT();
	// now derivatives
	double dd = tp.dDdP()[ipar]*dpar;
	double dt = tp.dTdP()[ipar]*dpar;
	dtpoca[ipar]->SetPoint(istep,xd-refd,dd);
	ttpoca[ipar]->SetPoint(istep,xt-reft,dt);
      }
    }
  }
  for(int ipar=0;ipar<lhel.npars_;ipar++){
    dtpcan->cd(ipar+1);
    dtpoca[ipar]->Draw("A*");
  }
  for(int ipar=0;ipar<lhel.npars_;ipar++){
    ttpcan->cd(ipar+1);
    ttpoca[ipar]->Draw("A*");
  }
  dtpcan->Write();
  ttpcan->Write();
  tpfile.Write();
  tpfile.Close();
  return 0;
}


