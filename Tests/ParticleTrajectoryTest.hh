// 
// test basic functions of the ParticleTrajectory, using the LoopHelix class
//
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/General/BFieldMap.hh"
#include "KinKal/General/PhysicalConstants.h"

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
using KinKal::Line;

void print_usage() {
  printf("Usage: ParticleTrajectoryTest --changedir i --delta f --tstep f --nsteps i \n");
}

template <class KTRAJ>
int ParticleTrajectoryTest(int argc, char **argv) {
  using PKTRAJ = ParticleTrajectory<KTRAJ>;
  using PTCA = PiecewiseClosestApproach<KTRAJ,Line>;
  double mom(105.0), cost(0.7), phi(0.5);
  unsigned npts(50);
  int icharge(-1);
  int idir(0), nsteps(10);
  double pmass(0.511), oz(100.0), ot(0.0);
  double tstart(0.0), tstep(1.0);
  double delta(0.01); // fractional change
  double hlen(200.0); // half-length of the wire
  double vprop(0.8), gap(2.0);

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
  MomBasis::Direction tdir =static_cast<MomBasis::Direction>(idir);
  cout << "Testing ParticleTrajectory with "
    << nsteps << " kinks in " << MomBasis::directionName(tdir) << " direction of size "
    << delta << endl;

  // create a helix
  VEC3 bnom(0.0,0.0,1.0);
  UniformBFieldMap BF(bnom); // 1 Tesla
  VEC4 origin(0.0,0.0,oz,ot);
  double sint = sqrt(1.0-cost*cost);
  MOM4 momv(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost,pmass);
  //initial piece
  double tend = tstart + tstep;
  TimeRange range(tstart,tend);
  KTRAJ lhel(origin,momv,icharge,bnom,range);
  // create initial piecewise helix from this
  PKTRAJ ptraj(lhel);
  // append pieces
  for(int istep=0;istep < nsteps; istep++){
// use derivatives of last piece to define new piece
    KTRAJ const& back = ptraj.pieces().back();
    double tcomp = back.range().end();
    DVEC pder = back.momDeriv(tcomp,tdir);
    // create modified helix
    DVEC dvec = back.params().parameters() + delta*pder;
    range = TimeRange(ptraj.range().end(),ptraj.range().end()+tstep);
    Parameters pdata(dvec,back.params().covariance());
    KTRAJ endhel(pdata,back);
    endhel.range() = range;
    // test
    VEC4 backpos = back.position4(tcomp);
    VEC4 endpos = endhel.position4(tcomp);
    cout << "back position " << backpos << endl
    << " end position " << endpos << endl;
    // append this
//    bool added(false);
    ptraj.append(endhel);
    // compare positions and momenta
    auto pold = back.position3(tcomp);
    auto mold = back.momentum4(tcomp);
    auto pnew = endhel.position3(tcomp);
    auto mnew = endhel.momentum4(tcomp);
    double dmom = sqrt((mnew.Vect()-mold.Vect()).Mag2());
    double gap = sqrt((pnew-pold).Mag2());
    cout << "Kink at time " << tcomp << " Position gap = " << gap << " Momentum change = " << dmom << endl;
  }
  // prepend pieces
  for(int istep=0;istep < nsteps; istep++){
    KTRAJ const& front = ptraj.pieces().front();
    double tcomp = front.range().begin();
    DVEC pder = front.momDeriv(tcomp,tdir);
    // create modified helix
    DVEC dvec = front.params().parameters() + delta*pder;
    range = TimeRange(ptraj.range().begin()-tstep,ptraj.range().begin());
    Parameters pdata(dvec,front.params().covariance());
    KTRAJ endhel(pdata,front);
    endhel.range() = range;
    // test
    auto frontpos = front.position4(tcomp);
    auto endpos = endhel.position4(tcomp);
    cout << "front position " << frontpos << endl
    << " end position " << endpos << endl;
    ptraj.prepend(endhel);
    // compare positions and momenta
    auto pold = front.position3(tcomp);
    auto mold = front.momentum4(tcomp);
    auto pnew = endhel.position3(tcomp);
    auto mnew = endhel.momentum4(tcomp);
    double dmom = sqrt((mnew.Vect()-mold.Vect()).Mag2());
    double gap = sqrt((pnew-pold).Mag2());
    cout << "Kink at time " << tcomp << " Position gap = " << gap << " Momentum change = " << dmom << endl;
  }
  double largest, average;
  size_t igap;
  ptraj.gaps(largest, igap, average);
  cout << "Final piece traj with " << ptraj.pieces().size() << " pieces and largest gap = "
  << largest << " average gap = " << average << endl;

// draw each piece of the piecetraj
  char fname[100];
  snprintf(fname,100,"ParticleTrajectory_%s_%2.2f.root",MomBasis::directionName(tdir).c_str(),delta);
  TFile pkfile((KTRAJ::trajName()+fname).c_str(),"RECREATE");
  TCanvas* pttcan = new TCanvas("pttcan","PieceKTRAJ",1000,1000);
  std::vector<TPolyLine3D*> plhel;
  int icolor(kBlue);
  for(auto const& piece : ptraj.pieces()) {
    plhel.push_back(new TPolyLine3D(npts));
    plhel.back()->SetLineColor(icolor);
    if(icolor == kBlue)
      icolor = kRed;
    else if(icolor == kRed)
      icolor = kBlue;
    double tstart = piece.range().begin();
    double ts = (piece.range().end()-piece.range().begin())/(npts-1);
    VEC3 ppos;
    for(unsigned ipt=0;ipt<npts;ipt++){
      double t = tstart + ipt*ts;
      ppos = piece.position3(t);
      plhel.back()->SetPoint(ipt,ppos.X(),ppos.Y(),ppos.Z());
    }
    plhel.back()->Draw();
  }
// now draw using the PTraj
  unsigned np = npts*ptraj.pieces().size();
  TPolyLine3D* all = new TPolyLine3D(np);
  all->SetLineColor(kYellow);
  all->SetLineStyle(kDotted);
  double ts = (ptraj.range().end()-ptraj.range().begin())/(np-1);
  for(unsigned ip=0;ip<np;ip++){
  double tp = ptraj.range().begin() + ip*ts;
    VEC3 ppos = ptraj.position3(tp);
      all->SetPoint(ip,ppos.X(),ppos.Y(),ppos.Z());
  }
  all->Draw();

  // draw the origin and axes
  TAxis3D* rulers = new TAxis3D();
  rulers->GetXaxis()->SetAxisColor(kBlue);
  rulers->GetXaxis()->SetLabelColor(kBlue);
  rulers->GetYaxis()->SetAxisColor(kCyan);
  rulers->GetYaxis()->SetLabelColor(kCyan);
  rulers->GetZaxis()->SetAxisColor(kOrange);
  rulers->GetZaxis()->SetLabelColor(kOrange);
  rulers->Draw();

// ClosestApproach test:
  VEC3 midpos, middir;
  midpos = ptraj.position3(ptraj.range().mid());
  middir = ptraj.direction(ptraj.range().mid());
  double lhphi = atan2(midpos.Y(),midpos.X());
  VEC3 rdir(cos(lhphi),sin(lhphi),0.0); // radial perpendicular to the helix
  double dphi = atan2(middir.Y(),middir.X());
  VEC3 pdir(-sin(dphi),cos(dphi),0.0);
  double pspeed = CLHEP::c_light*vprop; // vprop is relative to c
  VEC3 pvel = pdir*pspeed;
  TimeRange prange(ptraj.range().mid()-hlen/pspeed, ptraj.range().mid()+hlen/pspeed);
  // shift the position
  VEC3 lpos = midpos + gap*rdir;
  double wlen(1000.0);
  Line tline(lpos, ptraj.range().mid(), pvel, wlen);
  // create ClosestApproach from these
  CAHint tphint(ptraj.range().mid(),0.0);
  PTCA tp(ptraj,tline, tphint,1e-8);
  cout << "ClosestApproach status " << tp.statusName() << " doca " << tp.doca() << " dt " << tp.deltaT() << endl;
  VEC3 thpos, tlpos;
  thpos = tp.particlePoca().Vect();
  tlpos = tp.sensorPoca().Vect();
  double refd = tp.doca();
  cout << " Helix Pos " << midpos << " ClosestApproach KTRAJ pos " << thpos << " ClosestApproach Line pos " << tlpos << endl;
  cout << " ClosestApproach particlePoca " << tp.particlePoca() << " ClosestApproach sensorPoca " << tp.sensorPoca()  << " DOCA " << refd << endl;
  if(tp.status() == ClosestApproachData::converged) {
    // draw the line and ClosestApproach
    TPolyLine3D* line = new TPolyLine3D(2);
    VEC3 plow, phigh;
    plow = tline.position3(tline.range().begin());
    phigh = tline.position3(tline.range().end());
    line->SetPoint(0,plow.X(),plow.Y(), plow.Z());
    line->SetPoint(1,phigh.X(),phigh.Y(), phigh.Z());
    line->SetLineColor(kOrange);
    line->Draw();
    TPolyLine3D* poca = new TPolyLine3D(2);
    poca->SetPoint(0,tp.particlePoca().X() ,tp.particlePoca().Y() ,tp.particlePoca().Z());
    poca->SetPoint(1,tp.sensorPoca().X() ,tp.sensorPoca().Y() ,tp.sensorPoca().Z());
    poca->SetLineColor(kBlack);
    poca->Draw();
  }

  cout << "ClosestApproach dDdP" << tp.dDdP() << " dTdP " << tp.dTdP() << endl;
 
  pttcan->Write();

  pkfile.Write();
  pkfile.Close();
  // 
  return 0;
}


