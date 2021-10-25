//
// test basic functions of BFieldMap class
//
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Trajectory/ClosestApproach.hh"
#include "KinKal/General/BFieldMap.hh"
#include "KinKal/General/Vectors.hh"
#include "KinKal/General/PhysicalConstants.h"
#include "KinKal/Tests/ToyMC.hh"

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
  printf("Usage: BFieldMapTest  --momentum f --charge i --dBz f --dBx f --dBy f --Bgrad f --Tol f n");
}

template <class KTRAJ>
int BFieldMapTest(int argc, char **argv) {
  using PKTRAJ = ParticleTrajectory<KTRAJ>;
  using HIT = Hit<KTRAJ>;
  using HITPTR = std::shared_ptr<HIT>;
  using EXING = ElementXing<KTRAJ>;
  using EXINGPTR = std::shared_ptr<EXING>;
  using EXINGCOL = std::vector<EXINGPTR>;
  using HITCOL = std::vector<HITPTR>;
  double mom(105.0);
  int icharge(-1);
  int iseed(124223);
  double pmass(0.511);
  BFieldMap *BF(0);
  double Bgrad(0.0), dBx(0.0), dBy(0.0), dBz(0.0);
  double tol(0.1);
  double zrange(3000.0); // tracker dimension

  static struct option long_options[] = {
    {"momentum",     required_argument, 0, 'm' },
    {"charge",     required_argument, 0, 'q'  },
    {"dBx",     required_argument, 0, 'x'  },
    {"dBy",     required_argument, 0, 'y'  },
    {"dBz",     required_argument, 0, 'Z'  },
    {"Bgrad",     required_argument, 0, 'g'  },
    {"Tol",     required_argument, 0, 't'  },
    {NULL, 0,0,0}
  };
  int opt;
  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"",
	  long_options, &long_index )) != -1) {
    switch (opt) {
      case 'm' : mom = atof(optarg);
		 break;
      case 'x' : dBx = atof(optarg);
		 break;
      case 'y' : dBy = atof(optarg);
		 break;
      case 'Z' : dBz = atof(optarg);
		 break;
      case 'g' : Bgrad = atof(optarg);
		 break;
      case 't' : tol = atof(optarg);
		 break;
      case 'q' : icharge = atoi(optarg);
		 break;
      default: print_usage();
	       exit(EXIT_FAILURE);
    }
  }
  VEC3 bnom(0.0,0.0,1.0);
  VEC3 bsim;
  if(Bgrad != 0){
    BF = new GradientBFieldMap(1.0-0.5*Bgrad,1.0+0.5*Bgrad,-0.5*zrange,0.5*zrange);
    bnom = BF->fieldVect(VEC3(0.0,0.0,0.0));
  } else {
    VEC3 bsim(dBx,dBy,1.0+dBz);
    BF = new UniformBFieldMap(bsim);
    bnom = VEC3(0.0,0.0,1.0);
  }
    // first, create a traj based on the actual field at this point
  KKTest::ToyMC<KTRAJ> toy(*BF, mom, icharge, zrange, iseed, 0, false, false,false, -1.0, pmass );
  toy.setTolerance(tol);
  PKTRAJ tptraj;
  HITCOL thits;
  EXINGCOL dxings;
  toy.simulateParticle(tptraj, thits, dxings);
  // then, create a piecetraj around the nominal field with corrections,
  auto pos = tptraj.position4(tptraj.range().begin());
  auto momv = tptraj.momentum4(pos.T());
  KTRAJ start(pos,momv,icharge,bnom,tptraj.range());
  PKTRAJ xptraj(start);
  PKTRAJ lptraj(start);
  // see how far I can go within tolerance
  TimeRange prange = start.range();
  do {
    auto const& piece = xptraj.back();
    prange.end() = BF->rangeInTolerance(piece,prange.begin(), tol);
// integrate the momentum change over this range
    VEC3 dp = BF->integrate(piece,prange);
    // approximate change in position
//    VEC3 dpos = 0.5*dp*piece.speed(prange.mid())*prange.range()/piece.momentum4(prange.mid());
    // create a new trajectory piece at this point, correcting for the momentum change
    pos = piece.position4(prange.end());
    momv = piece.momentum4(pos.T());
//    cout << "BFieldMap integral dP " << dp.R() << " dpos " << dpos.R()  << " range " << prange.range() << " pos " << pos << endl;
    prange = TimeRange(prange.end(),std::max(prange.end()+double(0.1),xptraj.range().end()));
    // clumsy code to modify a vector
    auto mom = momv.Vect();
    mom += dp;
    momv.SetPx(mom.X()); momv.SetPy(mom.Y()); momv.SetPz(mom.Z());
    KTRAJ xnew(pos,momv,icharge,bnom,prange);
    // append this
    xptraj.append(xnew);
    double gap = xptraj.gap(xptraj.pieces().size()-1);
    if(gap > 1e-6) cout << "Apended traj gap " << gap << endl;
    // same thing, but using the linear parameter corrections
    VEC3 t1hat = lptraj.back().direction(pos.T(),MomBasis::perpdir_);
    VEC3 t2hat = lptraj.back().direction(pos.T(),MomBasis::phidir_);
    DVEC dpdt1 = lptraj.back().momDeriv(pos.T(),MomBasis::perpdir_);
    DVEC dpdt2 = lptraj.back().momDeriv(pos.T(),MomBasis::phidir_);
    auto dpfrac = dp/mom.R();
    DVEC dpars = dpfrac.Dot(t1hat)*dpdt1 + dpfrac.Dot(t2hat)*dpdt2;
    KTRAJ lnew = lptraj.back();
    lnew.params() += dpars;
    lnew.setRange(prange);
    lptraj.append(lnew);
    gap = xptraj.gap(xptraj.pieces().size()-1);
    if(gap > 1e-6) cout << "Apended traj gap " << gap << endl;
  }  while(prange.begin() < tptraj.range().end());
  // test integrating the field over the corrected trajectories: this should be small
  VEC3 tdp, xdp, ldp, ndp;
  tdp = BF->integrate( tptraj, tptraj.range());
  xdp = BF->integrate( xptraj, xptraj.range());
  ldp = BF->integrate( lptraj, lptraj.range());
  ndp = BF->integrate( start, start.range());
  cout << "TTraj " << tptraj << " integral " << tdp << endl;
  cout << "XTraj " << xptraj << " integral " << xdp << endl;
  cout << "LTraj " << lptraj << " integral " << ldp << endl;
  cout << "Nominal " << start << " integral " << ndp << endl;

// setup histograms
  TFile tpfile((KTRAJ::trajName()+"BFieldMap.root").c_str(),"RECREATE");
  TH1F *dxmomt1, *dxmomt2, *dxmommd;
  TH1F *dlmomt1, *dlmomt2, *dlmommd;
  TH1F *dxpost1, *dxpost2, *dxposmd;
  TH1F *dlpost1, *dlpost2, *dlposmd;
  dxmomt1 = new TH1F("dxmomt1","#Delta mom, T1 Dir",100,-2.0,2.0);
  dxmomt2 = new TH1F("dxmomt2","#Delta mom, T2 Dir",100,-2.0,2.0);
  dxmommd = new TH1F("dxmommd","#Delta mom, M Dir",100,-2.0,2.0);

  dlmomt1 = new TH1F("dlmomt1","#Delta mom, T1 Dir",100,-2.0,2.0);
  dlmomt2 = new TH1F("dlmomt2","#Delta mom, T2 Dir",100,-2.0,2.0);
  dlmommd = new TH1F("dlmommd","#Delta mom, M Dir",100,-2.0,2.0);

  dxpost1 = new TH1F("dxpost1","#Delta Pos, T1 Dir",100,-5.0,5.0);
  dxpost2 = new TH1F("dxpost2","#Delta Pos, T2 Dir",100,-5.0,5.0);
  dxposmd = new TH1F("dxposmd","#Delta Pos, M Dir",100,-5.0,5.0);

  dlpost1 = new TH1F("dlpost1","#Delta Pos, T1 Dir",100,-5.0,5.0);
  dlpost2 = new TH1F("dlpost2","#Delta Pos, T2 Dir",100,-5.0,5.0);
  dlposmd = new TH1F("dlposmd","#Delta Pos, M Dir",100,-5.0,5.0);

  // draw the true helix
  TCanvas* hcan = new TCanvas("hcan","Traj",1000,1000);
  TPolyLine3D* ttrue = new TPolyLine3D(200);
  TPolyLine3D* tpx = new TPolyLine3D(200);
  TPolyLine3D* tpl = new TPolyLine3D(200);
  TPolyLine3D* tnom = new TPolyLine3D(200);
  VEC3 tpos, xpos, lpos, npos;
  VEC3 tvel, xvel, lvel;
  VEC3 t1dir, t2dir, mdir;
  VEC3 tmom, xmom, lmom;
  double tstep = tptraj.range().range()/200.0;
  for(int istep=0;istep<201;++istep){
  // compute the position from the time
    double ttime = tptraj.range().begin() + tstep*istep;
    tpos = tptraj.position3(ttime);
    t1dir = tptraj.direction(ttime,MomBasis::perpdir_);
    t2dir = tptraj.direction(ttime,MomBasis::phidir_);
    mdir = tptraj.direction(ttime,MomBasis::momdir_);
    ttrue->SetPoint(istep, tpos.X(), tpos.Y(), tpos.Z());

    xpos = xptraj.position3(ttime);
    tpx->SetPoint(istep, xpos.X(), xpos.Y(), xpos.Z());

    lpos = lptraj.position3(ttime);
    tpl->SetPoint(istep, lpos.X(), lpos.Y(), lpos.Z());

    npos = start.position3(ttime);

    tnom->SetPoint(istep, npos.X(), npos.Y(), npos.Z());

    dxpost1->Fill( (xpos-tpos).Dot(t1dir));
    dxpost2->Fill( (xpos-tpos).Dot(t2dir));
    dxposmd->Fill( (xpos-tpos).Dot(mdir));

    dlpost1->Fill( (lpos-tpos).Dot(t1dir));
    dlpost2->Fill( (lpos-tpos).Dot(t2dir));
    dlposmd->Fill( (lpos-tpos).Dot(mdir));

    tmom = tptraj.momentum3(ttime);
    xmom = xptraj.momentum3(ttime);
    lmom = lptraj.momentum3(ttime);

    dxmomt1->Fill( (xmom-tmom).Dot(t1dir));
    dxmomt2->Fill( (xmom-tmom).Dot(t2dir));
    dxmommd->Fill( (xmom-tmom).Dot(mdir));

    dlmomt1->Fill( (lmom-tmom).Dot(t1dir));
    dlmomt2->Fill( (lmom-tmom).Dot(t2dir));
    dlmommd->Fill( (lmom-tmom).Dot(mdir));
  }
  // draw the true trajectory
  ttrue->SetLineColor(kBlue);
  ttrue->Draw();
  // draw the nominal (unadjusted) trajectory
  tnom->SetLineColor(kBlack);
  tnom->SetLineStyle(9);
  tnom->Draw();
  // draw the 'exact' (rotated) adjusted trajectory
  tpx->SetLineColor(kGreen);
  tpx->SetLineStyle(6);
  tpx->Draw();
  // linear adjusted trajectory
  tpl->SetLineColor(kRed);
  tpl->SetLineStyle(7);
  tpl->Draw();

  // draw the origin and axes
  TAxis3D* rulers = new TAxis3D();
  rulers->GetXaxis()->SetAxisColor(kBlue);
  rulers->GetXaxis()->SetLabelColor(kBlue);
  rulers->GetYaxis()->SetAxisColor(kCyan);
  rulers->GetYaxis()->SetLabelColor(kCyan);
  rulers->GetZaxis()->SetAxisColor(kOrange);
  rulers->GetZaxis()->SetLabelColor(kOrange);
  rulers->Draw();
  TLegend* hleg = new TLegend(0.7,0.7,1.0,1.0);
  hleg->AddEntry(ttrue,"True Trajectory","L");
  hleg->AddEntry(tnom,"Nominal Trajectory","L");
  hleg->AddEntry(tpx,"Rotated Correction Trajectory","L");
  hleg->AddEntry(tpl,"Fixed Correction Trajectory","L");
  hleg->Draw();
  hcan->Draw();
  hcan->Write();

  TCanvas* dxcan = new TCanvas("dxcan","dxcan",800,1200);
  dxcan->Divide(3,2);
  dxcan->cd(1);
  dxpost1->Draw();
  dxcan->cd(2);
  dxpost2->Draw();
  dxcan->cd(3);
  dxposmd->Draw();
  dxcan->cd(4);
  dxmomt1->Draw();
  dxcan->cd(5);
  dxmomt2->Draw();
  dxcan->cd(6);
  dxmommd->Draw();
  dxcan->Draw();
  dxcan->Write();

  TCanvas* dlcan = new TCanvas("dlcan","dlcan",800,1200);
  dlcan->Divide(3,2);
  dlcan->cd(1);
  dlpost1->Draw();
  dlcan->cd(2);
  dlpost2->Draw();
  dlcan->cd(3);
  dlposmd->Draw();
  dlcan->cd(4);
  dlmomt1->Draw();
  dlcan->cd(5);
  dlmomt2->Draw();
  dlcan->cd(6);
  dlmommd->Draw();
  dlcan->Draw();
  dlcan->Write();


// test the BFieldMap rotation by expressing a simple helix.  This doesn't simulate a true
// particle path, it just tests mechanics
// loop over the time ranges given by the 'sim' trajectory
  KTRAJ const& sktraj = tptraj.front();
  PKTRAJ rsktraj(sktraj);
  for (auto const& piece : tptraj.pieces()) {
    // rotate the parameters at the end of this piece to form the next.  Sample B in the middle of that range
    KTRAJ newpiece(piece,BF->fieldVect(sktraj.position3(piece.range().mid())),piece.range().begin());
    rsktraj.append(newpiece);
  }
  // draw the trajs
  TCanvas* hcandb = new TCanvas("hcandb","#Delta B Traj",1000,1000);
  TPolyLine3D* original = new TPolyLine3D(200);
  TPolyLine3D* rotated = new TPolyLine3D(200);
  tstep = tptraj.range().range()/200.0;
  for(int istep=0;istep<201;++istep){
    double tplot = sktraj.range().begin()+tstep*istep;
    auto opos = sktraj.position3(tplot);
    original->SetPoint(istep, opos.X(), opos.Y(), opos.Z());
    auto rpos = rsktraj.position3(tplot);
    rotated->SetPoint(istep, rpos.X(), rpos.Y(), rpos.Z());
  }
  original->SetLineColor(kBlack);
  original->SetLineStyle(9);
  original->Draw();
  rotated->SetLineColor(kRed);
  rotated->SetLineStyle(6);
  rotated->Draw();
  TLegend* bleg = new TLegend(0.7,0.7,1.0,1.0);
  bleg->AddEntry(original,"Original Trajectory","L");
  bleg->AddEntry(rotated,"Rotated Trajectory","L");
  bleg->Draw();
  hcandb->Draw();
  hcandb->Write();

  return 0;
}


