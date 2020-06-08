// 
// test basic functions of the PKTraj, using the LHelix class
//
#include "KinKal/PKTraj.hh"
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
  printf("Usage: PKTrajTest --changedir i --delta f --tstep f --nsteps i \n");
}

template <class KTRAJ>
int PKTrajTest(int argc, char **argv) {
  typedef PKTraj<KTRAJ> PKTRAJ;
  typedef typename KTRAJ::DVEC DVEC;
  typedef typename KTRAJ::PDATA PDATA; 
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
  LocalBasis::LocDir tdir =static_cast<LocalBasis::LocDir>(idir);
  cout << "Testing PKTraj with "
    << nsteps << " kinks in " << LocalBasis::directionName(tdir) << " direction of size "
    << delta << endl;

  // create a helix
  Vec3 bnom(0.0,0.0,1.0);
  UniformBField BF(bnom); // 1 Tesla
  Vec4 origin(0.0,0.0,oz,ot);
  double sint = sqrt(1.0-cost*cost);
  Mom4 momv(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost,pmass);
  //initial piece
  double tend = tstart + tstep;
  TRange range(tstart,tend);
  KTRAJ lhel(origin,momv,icharge,bnom,range);
  // create initial piecewise helix from this
  PKTRAJ ptraj(lhel);
  // append pieces
  for(int istep=0;istep < nsteps; istep++){
// use derivatives of last piece to define new piece
    KTRAJ const& back = ptraj.pieces().back();
    double tcomp = back.range().high();
    DVEC pder(back.momDeriv(tcomp,tdir));
    // create modified helix
    DVEC dvec1 = back.params().parameters();
    dvec1 += delta*pder;

    range = TRange(ptraj.range().high(),ptraj.range().high()+tstep);
    PDATA pdata(dvec1,back.params().covariance());
    KTRAJ endhel(pdata,back);
    endhel.range() = range;
    // test
    Vec4 backpos, endpos;
    backpos.SetE(tcomp);
    endpos.SetE(tcomp);
    back.position(backpos);
    endhel.position(endpos);
    cout << "back position " << backpos << endl
    << " end position " << endpos << endl;
    // append this
//    bool added(false);
    ptraj.append(endhel);
    // compare positions and momenta
    Vec3 pold, pnew;
    Mom4 mold, mnew;
    pold = back.position(tcomp);
    mold = back.momentum(tcomp);
    pnew = endhel.position(tcomp);
    mnew = endhel.momentum(tcomp);
    double dmom = sqrt((mnew.Vect()-mold.Vect()).Mag2());
    double gap = sqrt((pnew-pold).Mag2());
    cout << "Kink at time " << tcomp << " Position gap = " << gap << " Momentum change = " << dmom << endl;
  }
  // prepend pieces
  for(int istep=0;istep < nsteps; istep++){
    KTRAJ const& front = ptraj.pieces().front();
    double tcomp = front.range().low();
    DVEC pder = front.momDeriv(tcomp,tdir);
    // create modified helix
    DVEC dvec = front.params().parameters();
    dvec += delta*pder;

    range = TRange(ptraj.range().low()-tstep,ptraj.range().low());
    PDATA pdata(dvec,front.params().covariance());
    KTRAJ endhel(pdata,front);
    endhel.range() = range;
    // test
    Vec4 frontpos, endpos;
    frontpos.SetE(tcomp);
    endpos.SetE(tcomp);
    front.position(frontpos);
    endhel.position(endpos);
    cout << "front position " << frontpos << endl
    << " end position " << endpos << endl;
    ptraj.prepend(endhel);
    // compare positions and momenta
    Vec3 pold, pnew;
    Mom4 mold, mnew;
    pold = front.position(tcomp);
    mold = front.momentum(tcomp);
    pnew = endhel.position(tcomp);
    mnew = endhel.momentum(tcomp);
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
  snprintf(fname,100,"PKTraj_%s_%2.2f.root",LocalBasis::directionName(tdir).c_str(),delta);
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
    double tstart = piece.range().low();
    double ts = (piece.range().high()-piece.range().low())/(npts-1);
    Vec3 ppos;
    for(unsigned ipt=0;ipt<npts;ipt++){
      double t = tstart + ipt*ts;
      ppos = piece.position(t);
      plhel.back()->SetPoint(ipt,ppos.X(),ppos.Y(),ppos.Z());
    }
    plhel.back()->Draw();
  }
// now draw using the PTraj
  unsigned np = npts*ptraj.pieces().size();
  TPolyLine3D* all = new TPolyLine3D(np);
  all->SetLineColor(kYellow);
  all->SetLineStyle(kDotted);
  double ts = (ptraj.range().high()-ptraj.range().low())/(np-1);
  for(unsigned ip=0;ip<np;ip++){
  double tp = ptraj.range().low() + ip*ts;
    Vec3 ppos = ptraj.position(tp);
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

// TPoca test:
  Vec3 midpos, middir;
  midpos = ptraj.position(ptraj.range().mid());
  middir = ptraj.direction(ptraj.range().mid());
  double lhphi = atan2(midpos.Y(),midpos.X());
  Vec3 rdir(cos(lhphi),sin(lhphi),0.0); // radial perpendicular to the helix
  double dphi = atan2(middir.Y(),middir.X());
  Vec3 pdir(-sin(dphi),cos(dphi),0.0);
  double pspeed = CLHEP::c_light*vprop; // vprop is relative to c
  Vec3 pvel = pdir*pspeed;
  TRange prange(ptraj.range().mid()-hlen/pspeed, ptraj.range().mid()+hlen/pspeed);
  // shift the position
  Vec3 lpos = midpos + gap*rdir;
  TLine tline(lpos, pvel,ptraj.range().mid(),prange);
  // create TPoca from these
  TPoca<PKTRAJ,TLine> tp(ptraj,tline);
  cout << "TPoca status " << tp.statusName() << " doca " << tp.doca() << " dt " << tp.deltaT() << endl;
  Vec3 thpos, tlpos;
  thpos = tp.particleTraj().position(tp.particleToca());
  tlpos = tp.sensorTraj().position(tp.sensorToca());
  double refd = tp.doca();
  cout << " Helix Pos " << midpos << " TPoca KTRAJ pos " << thpos << " TPoca TLine pos " << tlpos << endl;
  cout << " TPoca particlePoca " << tp.particlePoca() << " TPoca sensorPoca " << tp.sensorPoca()  << " DOCA " << refd << endl;
  if(tp.status() == TPocaBase::converged) {
    // draw the line and TPoca
    TPolyLine3D* line = new TPolyLine3D(2);
    Vec3 plow, phigh;
    plow = tline.position(tline.range().low());
    phigh = tline.position(tline.range().high());
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

  // now derivatives
  TPoca<PKTRAJ,TLine> tdp(tp);
  cout << "TPoca dDdP" << tdp.dDdP() << " dTdP " << tdp.dTdP() << endl;
 
  pttcan->Write();

  pkfile.Write();
  pkfile.Close();
  // 
  return 0;
}


