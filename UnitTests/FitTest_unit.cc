// 
// ToyMC test of fitting an LHelix-based KKTrk
//
#include "KinKal/PKTraj.hh"
#include "KinKal/LHelix.hh"
#include "KinKal/TLine.hh"
#include "KinKal/TPoca.hh"
#include "KinKal/StrawHit.hh"
#include "KinKal/Context.hh"
#include "KinKal/Types.hh"
#include "KinKal/Constants.hh"
#include "KinKal/KKHit.hh"
#include "KinKal/KKTrk.hh"

#include <iostream>
#include <getopt.h>
#include <typeinfo>

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
typedef KinKal::KKTrk<LHelix> KKTRK;
// ugly global variables
double zrange(3000.0), rmax(800.0); // tracker dimension
double sprop(0.8*KinKal::c_), sdrift(0.06), rstraw(2.5);
double sigt(0.3); // drift time resolution in ns

void print_usage() {
  printf("Usage: FitTest  --momentum f --costheta f --azimuth f --particle i --charge i --zrange f --nhits i --hres f --seed i --escale f\n");
}

// helper function
KinKal::TLine GenerateStraw(PLHelix const& traj, double htime, TRandom* TR) {
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
//  cout << "Generating hit at position " << dpos << endl;
  // find the ends where this reaches the cylinder
  double drho = dpos.Rho();
  double ddot = dpos.Dot(sdir);
  double shlen = sqrt(rmax*rmax + ddot*ddot - drho*drho); // straw half-length
  // choose the closest end to be the measurement end
  double rprop = (fabs(-ddot-shlen) < fabs(-ddot+shlen)) ?  -ddot-shlen : -ddot+shlen;
  // sign propagation velocity away from the measurement
  if(rprop>0){
    sdir *= -1.0;
    rprop *= -1.0;
  }

  Vec3 mpos = dpos - sdir*rprop;
  Vec3 vprop = sdir*sprop;
  // measured time is after propagation and drift
  double tmeas = hpos.T() + fabs(rprop)/sprop + fabs(rdrift)/sdrift;
  // smear measurement time
  tmeas = TR->Gaus(tmeas,sigt);
  // measurement time is the longest time
  TRange trange(tmeas-2.0*shlen/sprop,tmeas);
  // construct the trajectory for this hit
  return TLine(mpos,vprop,tmeas,trange);
}


int main(int argc, char **argv) {
  int opt;
  double mom(105.0), cost(0.7), phi(0.5);
  double masses[5]={0.511,105.66,139.57, 493.68, 938.0};
  int imass(0), icharge(-1);
  double pmass;
  float escale(1.0);;
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
    {"hres",     required_argument, 0, 'h'  },
    {"nhits",     required_argument, 0, 'n'  },
    {"escale",     required_argument, 0, 'e'  },
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
      case 'h' : sigt = atof(optarg);
		 break;
      case 'n' : nhits = atoi(optarg);
		 break;
      case 's' : iseed = atoi(optarg);
		 break;
      case 'e' : escale = atof(optarg);
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
  CVD2T d2t(sdrift,sigt*sigt); 
  Vec4 origin(0.0,0.0,0.0,0.0);
  float sint = sqrt(1.0-cost*cost);
  Mom4 momv(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost,pmass);
  LHelix lhel(origin,momv,icharge,context);
  cout << "True " << lhel << endl; 
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
  std::vector<StrawHit> shits; // this owns the hits
  cout << "Creating " << nhits << " hits " << endl;
  for(size_t ihit=0; ihit<nhits; ihit++){
    double htime = plhel.range().low() + ihit*dt;
    auto tline = GenerateStraw(plhel,htime,TR);
    TPoca<PLHelix,TLine> tp(plhel,tline);
    WireHit::LRAmbig ambig = tp.doca() < 0 ? WireHit::left : WireHit::right;
    // construct the hit from this trajectory
    StrawHit sh(tline,context,d2t,rstraw,ambig);
    shits.push_back(sh);
  }
  std::vector<const THit*> thits; // this references them as THits
  for(auto& shit : shits) thits.push_back(&shit);
  cout << "vector of hit points " << thits.size() << endl;
//  for(auto thit : thits) {
//    cout << "Hit at time " << thit->sensorTraj().range().mid() << endl;
//  }
// create the fit seed by randomizing the parameters
  double sigmas[6] = { 0.5, 0.5, 0.5, 0.5, 0.005, 5.0};
  auto seedpar = lhel.params().parameters();
  ROOT::Math::SMatrix<double, 6,6, ROOT::Math::MatRepSym<double,6> > seedcov = ROOT::Math::SMatrixIdentity();
  for(unsigned ipar=0;ipar < 6; ipar++){
    seedpar[ipar] += TR->Gaus(0.0,sigmas[ipar]*escale);
    seedcov[ipar][ipar] *= sigmas[ipar]*sigmas[ipar];
  }
  cout << "Createing seed helix params " << seedpar <<" covariance " << endl << seedcov << endl;
  LHelix seedhel(seedpar,seedcov,lhel.mass(),lhel.charge(),context,lhel.range());
// Create the KKTrk from these hits
  KKTRK kktrk(seedhel,thits,1e6);
  auto& effs = kktrk.effects();
// fit the track
  kktrk.fit();
  cout << "KKTrk fit status = " << kktrk.status() << endl;
  for(auto const& eff : effs) {
    cout << "Eff at time " << eff->time() << " status " << eff->status(TDir::forwards)  << " " << eff->status(TDir::backwards) << endl;
    auto ihit = dynamic_cast<const KKHit<LHelix>*>(eff.get());
    if(ihit != 0){
      cout << "Hit status " << ihit->poca().status() << " doca " << ihit->poca().doca() << endl;
      cout << "Hit residual = " << ihit->resid() << endl;
    }
  }
  exit(EXIT_SUCCESS);
}
