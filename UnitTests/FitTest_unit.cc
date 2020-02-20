// 
// ToyMC test of an LHelix-based fit 
//
#include "KinKal/PKTraj.hh"
#include "KinKal/LHelix.hh"
#include "KinKal/TLine.hh"
#include "KinKal/TPOCA.hh"
#include "KinKal/StrawHit.hh"
#include "KinKal/Context.hh"
#include "KinKal/Constants.hh"

#include <iostream>
#include <stdio.h>
#include <iostream>
#include <getopt.h>

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
#include "TF1.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TProfile2D.h"

using namespace KinKal;
using namespace std;
// avoid confusion with root
using KinKal::TLine;
typedef KinKal::PKTraj<LHelix> PLHelix;

void print_usage() {
  printf("Usage: FitTest  --momentum f --costheta f --azimuth f --particle i --charge i --zrange f --nhits i --seed i\n");
}

// helper function
KinKal::StrawHit GenerateStrawHit(PLHelix const& traj, Context const& context, D2T const& d2t,double htime, double sprop, double rmax, double sigt, double rstraw,TRandom* TR) {
  // start with the true helix position at this time
  Vec4 hpos; hpos.SetE(htime);
  traj.position(hpos);
  Vec3 hdir; traj.direction(htime,hdir);
  // generate a random direction for the straw 
  double eta = TR->Uniform(-M_PI,M_PI);
  Vec3 sdir(cos(eta),sin(eta),0.0);
  // generate a random drift perp to this and the trajectory
  double rdrift = TR->Uniform(-rstraw,rstraw);
  Vec3 drift = (sdir.Cross(hdir)).Unit();
  Vec3 dpos = hpos.Vect() + rdrift*drift;
  // find the ends where this reaches the cyldinder
  double drho = dpos.Rho();
  double ddot = dpos.Dot(sdir);
  double shlen = sqrt(rmax*rmax + ddot*ddot - drho*drho); // straw half-length
  // choose the closest end to be the measurement end
  double rprop;
  if(fabs(-ddot-shlen) < fabs(-ddot+shlen))
    rprop = -ddot-shlen;
  else
    rprop = -ddot+shlen;
  // sign propagation velocity away from the measurement
  if(rprop>0) sdir *= -1.0;
  Vec3 dd = dpos-hpos.Vect();
  Vec3 cross = sdir.Cross(hdir);
  WireHit::LRAmbig ambig = dd.Dot(cross) > 0 ? WireHit::left : WireHit::right;
  Vec3 mpos = dpos - sdir*rprop;
  Vec3 vprop = sdir*sprop;
  double tmeas = hpos.T() + rprop/sprop + rdrift/d2t.averageDriftSpeed();
  // smear measurement time
  tmeas = TR->Gaus(tmeas,sigt);
  // measurement time is the longest time
  TRange trange(tmeas-2.0*shlen/sprop,tmeas);
  // construct the trajectory for this hit
  TLine sline(mpos,vprop,tmeas,trange);
  // construct the hit from this trajectory
  return StrawHit(sline,context,d2t,rstraw,ambig);
}


int main(int argc, char **argv) {
  int opt;
  double mom(105.0), cost(0.7), phi(0.5);
  double masses[5]={0.511,105.66,139.57, 493.68, 938.0};
  int imass(0), icharge(-1);
  double pmass;
  double zrange(3000.0), rmax(800.0); // tracker dimension
  double sprop(0.8*KinKal::c_), sdrift(0.06), rstraw(2.5);
  double sigt(0.3); // drift time resolution in ns
  unsigned nhits(40);
  int iseed(124223);

  static struct option long_options[] = {
    {"momentum",     required_argument, 0, 'm' },
    {"costheta",     required_argument, 0, 'c'  },
    {"azimuth",     required_argument, 0, 'a'  },
    {"particle",     required_argument, 0, 'p'  },
    {"charge",     required_argument, 0, 'q'  },
    {"zrange",     required_argument, 0, 'z'  },
    {"seed",     required_argument, 0, 's'  },
    {"nhits",     required_argument, 0, 'n'  },
  };

  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"", 
	  long_options, &long_index )) != -1) {
    switch (opt) {
      case 'm' : mom = atof(optarg); 
		 break;
      case 'c' : cost = atof(optarg);
		 break;
      case 'a' : phi = atof(optarg);
		 break;
      case 'p' : imass = atoi(optarg);
		 break;
      case 'q' : icharge = atoi(optarg);
		 break;
      case 'z' : zrange = atof(optarg);
		 break;
      case 'n' : nhits = atoi(optarg);
		 break;
      case 's' : iseed = atoi(optarg);
		 break;
      default: print_usage(); 
	       exit(EXIT_FAILURE);
    }
  }

  TRandom* TR = new TRandom3(iseed);

  pmass = masses[imass];
// define the context
  UniformBField BF(1.0); // 1 Tesla
  Context context(BF);
  CVD2T d2t(sdrift,3.5); 
  Vec4 origin(0.0,0.0,0.0,0.0);
  float sint = sqrt(1.0-cost*cost);
  Mom4 momv(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost,pmass);
  LHelix lhel(origin,momv,icharge,context);
  PLHelix plhel(lhel);
  // truncate the range according to the Z range
  Vec3 vel; plhel.velocity(0.0,vel);
  plhel.range() = TRange(-0.5*zrange/vel.Z(),0.5*zrange/vel.Z());
  // generate material effects
  // FIXME!
  // divide time range
  double dt = plhel.range().range()/(nhits-1);
  cout << "True " << plhel << endl;
  // generate hits
  std::vector<StrawHit> hits;
  for(size_t ihit=0; ihit<hits.size(); ihit++){
    double htime = plhel.range().low() + ihit*dt;
    hits.push_back(GenerateStrawHit(plhel,context,d2t,htime,sprop,rmax,sigt,rstraw,TR));
  }
  exit(EXIT_SUCCESS);
}



