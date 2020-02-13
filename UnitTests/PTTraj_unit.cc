// 
// test basic functions of the PTTraj, using the LHelix class
//
#include "BTrk/KinKal/PKTraj.hh"
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
  unsigned npts(50);
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
  PKTraj<LHelix> ptraj(lhel);
  // append pieces
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
    range = TRange(ptraj.range().high(),ptraj.range().high()+tstep);
    LHelix endhel(dvec,back.params().mat(),back.mass(),back.charge(),context,range);
    // test
    Vec4 backpos, endpos;
    backpos.SetE(tcomp);
    endpos.SetE(tcomp);
    back.position(backpos);
    endhel.position(endpos);
    cout << "back position " << backpos << endl
    << " end position " << endpos << endl;
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
  // prepend pieces
  for(int istep=0;istep < nsteps; istep++){
    LHelix::PDer pder;
    LHelix const& front = ptraj.pieces().front();
    double tcomp = front.range().low();
    front.momDeriv(tdir,tcomp,pder);
    // create modified helix
    LHelix::TDATA::DVec dvec = front.params().vec();
    for(size_t ipar=0;ipar<6;ipar++)
      dvec[ipar] -= delta*pder[ipar][0];
    range = TRange(ptraj.range().low()-tstep,ptraj.range().low());
    LHelix endhel(dvec,front.params().mat(),front.mass(),front.charge(),context,range);
    // test
    Vec4 frontpos, endpos;
    frontpos.SetE(tcomp);
    endpos.SetE(tcomp);
    front.position(frontpos);
    endhel.position(endpos);
    cout << "front position " << frontpos << endl
    << " end position " << endpos << endl;
    bool added = ptraj.prepend(endhel);
    // compare positions and momenta
    Vec3 pold, pnew;
    Mom4 mold, mnew;
    front.position(tcomp,pold);
    front.momentum(tcomp,mold);
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

// draw each piece of the piecetraj
  TCanvas* pttcan = new TCanvas("pttcan","PieceLHelix",1000,1000);
  std::vector<TPolyLine3D*> plhel;
  int icolor(kBlue);
  for(auto const& piece : ptraj.pieces()) {
    plhel.push_back(new TPolyLine3D(npts));
    plhel.back()->SetLineColor(icolor);
    if(icolor == kBlue)
      icolor = kRed;
    else if(icolor == kRed)
      icolor = kBlue;
    double tstart = piece.range().low();
    double tstep = (piece.range().high()-piece.range().low())/(npts-1);
    Vec3 ppos;
    for(unsigned ipt=0;ipt<npts;ipt++){
      double t = tstart + ipt*tstep;
      piece.position(t,ppos);
      plhel.back()->SetPoint(ipt,ppos.X(),ppos.Y(),ppos.Z());
    }
    plhel.back()->Draw();
  }
  // draw the origin and axes
  TAxis3D* rulers = new TAxis3D();
  rulers->GetXaxis()->SetAxisColor(kBlue);
  rulers->GetXaxis()->SetLabelColor(kBlue);
  rulers->GetYaxis()->SetAxisColor(kCyan);
  rulers->GetYaxis()->SetLabelColor(kCyan);
  rulers->GetZaxis()->SetAxisColor(kOrange);
  rulers->GetZaxis()->SetLabelColor(kOrange);
  rulers->Draw();
  char fname[100];
  snprintf(fname,100,"PTTraj_%s_%2.2f.root",KTraj::directionName(tdir).c_str(),delta);
  pttcan->SaveAs(fname); 

  // 
  return 0;
}


