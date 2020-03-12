// 
// test basic functions of the PKTraj, using the LHelix class
//
#include "KinKal/PKTraj.hh"
#include "KinKal/LHelix.hh"
#include "KinKal/TLine.hh"
#include "KinKal/TPoca.hh"
#include "KinKal/Context.hh"
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
#include "TF1.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TProfile2D.h"

using namespace KinKal;
using namespace std;
// avoid confusion with root
using KinKal::TLine;

void print_usage() {
  printf("Usage: PKTraj --changedir i --delta f --tstep f --nsteps i \n");
}

int main(int argc, char **argv) {
  typedef PKTraj<LHelix> PLHelix;

  double mom(105.0), cost(0.7), phi(0.5);
  unsigned npts(50);
  int icharge(-1), idir(0), nsteps(10);
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
  KInter::MDir tdir =static_cast<KInter::MDir>(idir);
  cout << "Testing PKTraj with "
    << nsteps << " kinks in " << KInter::directionName(tdir) << " direction of size "
    << delta << endl;

  // create a helix
  UniformBField BF(1.0); // 1 Tesla
  Context context(BF);
  Vec4 origin(0.0,0.0,oz,ot);
  float sint = sqrt(1.0-cost*cost);
  Mom4 momv(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost,pmass);
  //initial piece
  double tend = tstart + tstep;
  TRange range(tstart,tend);
  LHelix lhel(origin,momv,icharge,context,range);
  // create initial piecewise helix from this
  PLHelix ptraj(lhel);
  // append pieces
  for(int istep=0;istep < nsteps; istep++){
// use derivatives of last piece to define new piece
    LHelix::PDER pder;
    LHelix const& back = ptraj.pieces().back();
    double tcomp = back.range().high();
    back.momDeriv(tdir,tcomp,pder);
    // create modified helix
    auto dvec = back.params().parameters() + delta*pder;
    range = TRange(ptraj.range().high(),ptraj.range().high()+tstep);
    LHelix endhel(dvec,back.params().covariance(),back.mass(),back.charge(),context,range);
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
    LHelix::PDER pder;
    LHelix const& front = ptraj.pieces().front();
    double tcomp = front.range().low();
    front.momDeriv(tdir,tcomp,pder);
    // create modified helix
    auto dvec = front.params().parameters() + delta*pder;
    range = TRange(ptraj.range().low()-tstep,ptraj.range().low());
    LHelix endhel(dvec,front.params().covariance(),front.mass(),front.charge(),context,range);
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
  char fname[100];
  snprintf(fname,100,"PKTraj_%s_%2.2f.root",KInter::directionName(tdir).c_str(),delta);
  TFile pkfile(fname,"RECREATE");
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
    double ts = (piece.range().high()-piece.range().low())/(npts-1);
    Vec3 ppos;
    for(unsigned ipt=0;ipt<npts;ipt++){
      double t = tstart + ipt*ts;
      piece.position(t,ppos);
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
    Vec3 ppos;
      ptraj.position(tp,ppos);
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
  ptraj.position(ptraj.range().mid(),midpos);
  ptraj.direction(ptraj.range().mid(),middir);
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
  TPoca<PLHelix,TLine> tp(ptraj,tline);
  cout << "TPoca status " << tp.statusName() << " doca " << tp.doca() << " dt " << tp.deltaT() << endl;
  Vec3 thpos, tlpos;
  tp.ttraj0().position(tp.poca0().T(),thpos);
  tp.ttraj1().position(tp.poca1().T(),tlpos);
  double refd = tp.doca();
  cout << " Helix Pos " << midpos << " TPoca LHelix pos " << thpos << " TPoca TLine pos " << tlpos << endl;
  cout << " TPoca poca0 " << tp.poca0() << " TPoca poca1 " << tp.poca1()  << " DOCA " << refd << endl;
  if(tp.status() == TPocaBase::converged) {
    // draw the line and TPoca
    TPolyLine3D* line = new TPolyLine3D(2);
    Vec3 plow, phigh;
    tline.position(tline.range().low(),plow);
    tline.position(tline.range().high(),phigh);
    line->SetPoint(0,plow.X(),plow.Y(), plow.Z());
    line->SetPoint(1,phigh.X(),phigh.Y(), phigh.Z());
    line->SetLineColor(kOrange);
    line->Draw();
    TPolyLine3D* poca = new TPolyLine3D(2);
    poca->SetPoint(0,tp.poca0().X() ,tp.poca0().Y() ,tp.poca0().Z());
    poca->SetPoint(1,tp.poca1().X() ,tp.poca1().Y() ,tp.poca1().Z());
    poca->SetLineColor(kBlack);
    poca->Draw();
  }

  // now derivatives
  TDPoca<PLHelix,TLine> tdp(tp);
  cout << "TDPoca dDdP" << tdp.dDdP() << " dTdP " << tdp.dTdP() << endl;
 
  pttcan->Write();

  pkfile.Write();
  pkfile.Close();
  // 
  return 0;
}


