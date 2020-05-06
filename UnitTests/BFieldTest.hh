// 
// test basic functions of BField class
//
#include "MatEnv/MatDBInfo.hh"
#include "MatEnv/DetMaterial.hh"
#include "KinKal/PKTraj.hh"
#include "KinKal/TLine.hh"
#include "KinKal/TPoca.hh"
#include "KinKal/StrawHit.hh"
#include "KinKal/StrawMat.hh"
#include "KinKal/BField.hh"
#include "KinKal/Vectors.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "UnitTests/ToyMC.hh"

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
  printf("Usage: BFieldTest  --momentum f --charge i --dBz f --dBx f --dBy f --Bgrad f --Tol f n");
}

template <class KTRAJ>
int BFieldTest(int argc, char **argv) {
  typedef KinKal::PKTraj<KTRAJ> PKTRAJ;
  typedef THit<KTRAJ> THIT;
  typedef std::shared_ptr<THIT> THITPTR;
  typedef DXing<KTRAJ> DXING;
  typedef std::shared_ptr<DXING> DXINGPTR;
  typedef std::vector<THITPTR> THITCOL;
  typedef vector<DXINGPTR> DXINGCOL;
  typedef typename KTRAJ::PDATA::DVEC DVEC;
  double mom(105.0);
  int icharge(-1);
  int iseed(124223);
  double pmass(0.511);
  BField *BF(0);
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
  Vec3 bnom(0.0,0.0,1.0);
  Vec3 bsim;
  if(Bgrad != 0){
    BF = new GradBField(1.0-0.5*Bgrad,1.0+0.5*Bgrad,-0.5*zrange,0.5*zrange);
    bnom = BF->fieldVect(Vec3(0.0,0.0,0.0));
  } else {
    Vec3 bsim(dBx,dBy,1.0+dBz);
    BF = new UniformBField(bsim);
    bnom = Vec3(0.0,0.0,1.0);
  }
    // first, create a traj based on the actual field at this point
  KKTest::ToyMC<KTRAJ> toy(*BF, mom, icharge, zrange, iseed, 0, false, false, -1.0, pmass );
  PKTRAJ tptraj;
  THITCOL thits;
  DXINGCOL dxings;
  toy.simulateParticle(tptraj, thits, dxings);
  // then, create a piecetraj around the nominal field with corrections,
  Vec4 pos; pos.SetE(tptraj.range().low());
  tptraj.position(pos);
  Mom4 momv = tptraj.momentum(pos.T());
  KTRAJ start(pos,momv,icharge,bnom,tptraj.range());
  PKTRAJ xptraj(start);
  PKTRAJ lptraj(start);
  // see how far I can go within tolerance
  TRange prange = start.range();
  do {
    auto const& piece = xptraj.back();
    piece.rangeInTolerance(prange,*BF, tol);
// integrate the momentum change over this range
    Vec3 dp;
    BF->integrate(piece,prange,dp);
    // approximate change in position
//    Vec3 dpos = 0.5*dp*piece.speed(prange.mid())*prange.range()/piece.momentum(prange.mid());
    // create a new trajectory piece at this point, correcting for the momentum change
    pos.SetE(prange.high());
    piece.position(pos);
    momv = piece.momentum(pos.T());
//    cout << "BField integral dP " << dp.R() << " dpos " << dpos.R()  << " range " << prange.range() << " pos " << pos << endl;
    prange = TRange(prange.high(),std::max(prange.high()+double(0.1),xptraj.range().high()));
    // clumsy code to modify a vector
    Vec3 mom = momv.Vect();
    mom += dp;
    momv.SetPx(mom.X()); momv.SetPy(mom.Y()); momv.SetPz(mom.Z());
    KTRAJ xnew(pos,momv,icharge,bnom,prange);
    // append this
    xptraj.append(xnew);
    double gap = xptraj.gap(xptraj.pieces().size()-1);
    if(gap > 1e-6) cout << "Apended traj gap " << gap << endl;
    // same thing, but using the linear parameter corrections
    Vec3 t1hat, t2hat;
    DVEC dpdt1, dpdt2;
    lptraj.back().momDeriv(pos.T(),LocalBasis::perpdir,dpdt1,t1hat);
    lptraj.back().momDeriv(pos.T(),LocalBasis::phidir,dpdt2,t2hat);
    auto dpfrac = dp/mom.R();
    DVEC dpars = dpfrac.Dot(t1hat)*dpdt1 + dpfrac.Dot(t2hat)*dpdt2;
    KTRAJ lnew = lptraj.back();
    lnew.params() += dpars;
    lnew.setRange(prange);
    lptraj.append(lnew);
    gap = xptraj.gap(xptraj.pieces().size()-1);
    if(gap > 1e-6) cout << "Apended traj gap " << gap << endl;
  }  while(prange.low() < tptraj.range().high());
  // test integrating the field over the corrected trajectories: this should be small
  Vec3 tdp, xdp, ldp, ndp;
  BF->integrate(tptraj, tptraj.range(),tdp);
  BF->integrate(xptraj, xptraj.range(),xdp);
  BF->integrate(lptraj, lptraj.range(),ldp);
  BF->integrate(start, start.range(),ndp);
  cout << "TTraj " << tptraj << " integral " << tdp << endl;
  cout << "XTraj " << xptraj << " integral " << xdp << endl;
  cout << "LTraj " << lptraj << " integral " << ldp << endl;
  cout << "Nominal " << start << " integral " << ndp << endl;

// setup histograms
  TFile tpfile("BField.root","RECREATE");
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
  Vec4 tpos, xpos, lpos, npos;
  Vec3 tvel, xvel, lvel;
  Vec3 t1dir, t2dir, mdir;
  Mom4 tmom, xmom, lmom;
  double tstep = tptraj.range().range()/200.0;
  for(int istep=0;istep<201;++istep){
  // compute the position from the time
    tpos.SetE(tptraj.range().low() + tstep*istep);
    tptraj.position(tpos);
    t1dir = tptraj.direction(tpos.T(),LocalBasis::perpdir);
    t2dir = tptraj.direction(tpos.T(),LocalBasis::phidir);
    mdir = tptraj.direction(tpos.T(),LocalBasis::momdir);
 
    ttrue->SetPoint(istep, tpos.X(), tpos.Y(), tpos.Z());
    xpos.SetE(tptraj.range().low() + tstep*istep);
    xptraj.position(xpos);
    tpx->SetPoint(istep, xpos.X(), xpos.Y(), xpos.Z());
    lpos.SetE(tptraj.range().low() + tstep*istep);
    lptraj.position(lpos);
    tpl->SetPoint(istep, lpos.X(), lpos.Y(), lpos.Z());
    npos.SetE(tptraj.range().low() + tstep*istep);
    start.position(npos);
    tnom->SetPoint(istep, npos.X(), npos.Y(), npos.Z());

    dxpost1->Fill( (xpos-tpos).Vect().Dot(t1dir));
    dxpost2->Fill( (xpos-tpos).Vect().Dot(t2dir));
    dxposmd->Fill( (xpos-tpos).Vect().Dot(mdir));

    dlpost1->Fill( (lpos-tpos).Vect().Dot(t1dir));
    dlpost2->Fill( (lpos-tpos).Vect().Dot(t2dir));
    dlposmd->Fill( (lpos-tpos).Vect().Dot(mdir));

    tmom = tptraj.momentum(tpos.T());
    xmom = xptraj.momentum(xpos.T());
    lmom = lptraj.momentum(lpos.T());
    dxmomt1->Fill( (xmom.Vect()-tmom.Vect()).Dot(t1dir));
    dxmomt2->Fill( (xmom.Vect()-tmom.Vect()).Dot(t2dir));
    dxmommd->Fill( (xmom.Vect()-tmom.Vect()).Dot(mdir));

    dlmomt1->Fill( (lmom.Vect()-tmom.Vect()).Dot(t1dir));
    dlmomt2->Fill( (lmom.Vect()-tmom.Vect()).Dot(t2dir));
    dlmommd->Fill( (lmom.Vect()-tmom.Vect()).Dot(mdir));
  }
  // draw the true helix
  ttrue->SetLineColor(kBlue);
  ttrue->Draw();
  // draw the nominal (unadjusted) helix
  tnom->SetLineColor(kBlack);
  tnom->SetLineStyle(9);
  tnom->Draw();
  // draw the adjusted helix
  tpx->SetLineColor(kGreen);
  tpx->SetLineStyle(6);
  tpx->Draw();
  //
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
  hleg->AddEntry(tpx,"Exact Correction Trajectory","L");
  hleg->AddEntry(tpl,"Linear Correction Trajectory","L");
  hleg->Draw();
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
 
  return 0;
}


