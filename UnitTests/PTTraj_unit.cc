// 
// test basic functions of the PTTraj, using the LHelix class
//
#include "BTrk/KinKal/PTTraj.hh"
#include "BTrk/KinKal/LHelix.hh"
#include "BTrk/KinKal/Context.hh"

#include <iostream>
#include <stdio.h>
#include <iostream>
#include <getopt.h>

#include "TH1F.h"
#include "TSystem.h"
#include "THelix.h"
#include "TPolyLine3D.h"
#include "TArrow.h"
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

void print_usage() {
  printf("Usage: PTTraj --changedir i --delta f --tstep f --nsteps i \n");
}

int main(int argc, char **argv) {
  double mom(105.0), cost(0.7), phi(0.5);
  int icharge(-1), idir(0), nsteps(10);
  double pmass(0.511), oz(100.0), ot(0.0);
  double tstart(0.0), tstep(1.0);
  double delta(0.01); // fractional change

  static struct option long_options[] = {
    {"changedir",     required_argument, 0, 'c'  },
    {"delta",     required_argument, 0, 'd'  },
    {"tstep",     required_argument, 0, 't'  },
    {"nsteps",     required_argument, 0, 'n'  }
  };

  int long_index =0;
  int opt;
  while ((opt = getopt_long_only(argc, argv,"", 
	  long_options, &long_index )) != -1) {
    switch (opt) {
      case 'c' : 
	idir = atoi(optarg);
	break;
      case 'd' :
	delta = atof(optarg);
	break;
      case 't' :
	tstep = atof(optarg);
	break;
      case 'n' :
	nsteps = atoi(optarg);
	break;
      default:
	print_usage(); 
	exit(EXIT_FAILURE);
    }
  }
  KTraj::trajdir tdir =static_cast<KTraj::trajdir>(idir);
  cout << "Testing PTTraj with "
    << nsteps << " kinks in " << KTraj::directionName(tdir) << " direction of size "
    << delta << endl;

  // create a helix
  Context context;
  context.Bz_ = 1.0; // 1 Tesla
  Vec4 origin(0.0,0.0,oz,ot);
  float sint = sqrt(1.0-cost*cost);
  Mom4 momv(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost,pmass);
  //initial piece
  double tend = tstart + tstep;
  TRange range(tstart,tend);
  LHelix lhel(origin,momv,icharge,context,range);
  // create piece
  PTTraj<LHelix> ptraj(lhel);
  // add pieces
  for(int istep=0;istep < nsteps; istep++){
// use derivatives of last piece to define new piece
    LHelix::PDer pder;
    LHelix const& back = ptraj.pieces().back();
    double tcomp = back.range().high();
    back.momDeriv(tdir,tcomp,pder);
    // create modified helix
    LHelix::TDATA::DVec dvec = back.params().vec();
    for(size_t ipar=0;ipar<6;ipar++)
      dvec[ipar] += delta*pder[ipar][0];
    range = TRange(range.high(),range.high()+tstep);
    LHelix endhel(dvec,back.params().mat(),back.mass(),back.charge(),context,range);
    // append this
    bool added = ptraj.append(endhel);
    // compare positions and momenta
    Vec3 pold, pnew;
    Mom4 mold, mnew;
    back.position(tcomp,pold);
    back.momentum(tcomp,mold);
    endhel.position(tcomp,pnew);
    endhel.momentum(tcomp,mnew);
    double dmom = sqrt((mnew.Vect()-mold.Vect()).Mag2());
    double gap = sqrt((pnew-pold).Mag2());
    cout << "Kink at time " << tcomp << " Position gap = " << gap << " Momentum change = " << dmom << endl;

    if(added)
      cout << "succeded to append LHelix" << endl;
    else
      cout << "failed to append LHelix" << endl;
  }
  double largest, average;
  size_t igap;
  ptraj.gaps(largest, igap, average);
  cout << "Final piece traj with " << ptraj.pieces().size() << " pieces and largest gap = "
  << largest << " average gap = " << average << endl;

  return 0;
}


