// 
// test basic functions of BField class
//
#include "MatEnv/MatDBInfo.hh"
#include "MatEnv/DetMaterial.hh"
#include "KinKal/PKTraj.hh"
#include "KinKal/LHelix.hh"
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
  printf("Usage: BField  --momentum f --charge i --dBz f --dBx f --dBy f --Bgrad f  \n");
}

int main(int argc, char **argv) {
  typedef LHelix KTRAJ;
  typedef KinKal::PKTraj<KTRAJ> PKTRAJ;
  typedef THit<KTRAJ> THIT;
  typedef std::shared_ptr<THIT> THITPTR;
  typedef DXing<KTRAJ> DXING;
  typedef std::shared_ptr<DXING> DXINGPTR;
  typedef std::vector<THITPTR> THITCOL;
  typedef vector<DXINGPTR> DXINGCOL;
  typedef typename KTRAJ::PDATA::DVEC DVEC;
  double mom(105.0);
  float tmin(-10.0), tmax(10.0);
  int icharge(-1);
  int iseed(124223);
  double pmass(0.511);
  BField *BF(0);
  float Bgrad(0.0), dBx(0.0), dBy(0.0), dBz(0.0);
  double zrange(3000.0); // tracker dimension

  static struct option long_options[] = {
    {"momentum",     required_argument, 0, 'm' },
    {"charge",     required_argument, 0, 'q'  },
    {"dBx",     required_argument, 0, 'x'  },
    {"dBy",     required_argument, 0, 'y'  },
    {"dBz",     required_argument, 0, 'Z'  },
    {"Bgrad",     required_argument, 0, 'g'  },
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
    BF->fieldVect(Vec3(0.0,0.0,0.0),bnom);
  } else {
    Vec3 bsim(dBx,dBy,1.0+dBz);
    BF = new UniformBField(bsim);
    bnom = Vec3(0.0,0.0,1.0);
  }
    // first, create a traj based on the actual field at this point
  KKTest::ToyMC<KTRAJ> toy(*BF, mom, icharge, zrange, iseed, 40, false, false, -1.0, pmass );
  PKTRAJ tptraj(TRange(tmin,tmax),pmass,icharge);
  THITCOL thits;
  DXINGCOL dxings;
  toy.simulateParticle(tptraj, thits, dxings);
  // then, create a piecetraj around the nominal field with corrections,
  Mom4 momv;
  Vec4 pos; pos.SetE(tptraj.range().low());
  tptraj.position(pos);
  tptraj.momentum(pos.T(),momv);
  KTRAJ start(pos,momv,icharge,bnom,tptraj.range());
  PKTRAJ xptraj(start);
  PKTRAJ lptraj(start);
  // see how far I can go within tolerance
  TRange prange = start.range();
  do {
    auto const& piece = xptraj.back();
    piece.rangeInTolerance(prange,*BF,0.1,0.1);
// integrate the momentum change over this range
    Vec3 dp;
    BF->integrate(piece,prange,dp);
    // create a new trajectory piece at this point, correcting for the momentum change
    pos.SetE(prange.high());
    piece.position(pos);
    piece.momentum(pos.T(),momv);
//    cout << "BField integral dP " << dp << " range " << prange.range() << " pos " << pos << endl;
    prange = TRange(prange.high(),std::max(prange.high()+float(0.1),xptraj.range().high()));
    // clumsy code to modify a vector
    Vec3 mom = momv.Vect();
    mom += dp;
    momv.SetPx(mom.X()); momv.SetPy(mom.Y()); momv.SetPz(mom.Z());
    KTRAJ xnew(pos,momv,icharge,bnom,prange);
    // append this
    xptraj.append(xnew);
    float gap = xptraj.gap(xptraj.pieces().size()-1);
    if(gap > 1e-6) cout << "Apended traj gap " << gap << endl;
    // same thing, but using the linear parameter corrections
    Vec3 t1hat, t2hat;
    DVEC dpdt1, dpdt2;
    lptraj.back().momDeriv(KInter::theta1,pos.T(),dpdt1,t1hat);
    lptraj.back().momDeriv(KInter::theta2,pos.T(),dpdt2,t2hat);
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
  TH1F* dxpos, *dxvel, *dxspeed, *dxmom;
  TH1F* dlpos, *dlvel, *dlspeed, *dlmom;
  dxpos = new TH1F("dxpos","Position difference",100,0.0,5.0);
  dxvel = new TH1F("dxvel","Velocity difference",100,0.0,10.0);
  dxspeed = new TH1F("dxspeed","Speed difference",100,0.0,1.0);
  dxmom = new TH1F("dxmom","momentum difference",100,0.0,2.0);

  dlpos = new TH1F("dlpos","Position difference",100,0.0,5.0);
  dlvel = new TH1F("dlvel","Velocity difference",100,0.0,10.0);
  dlspeed = new TH1F("dlspeed","Speed difference",100,0.0,1.0);
  dlmom = new TH1F("dlmom","momentum difference",100,0.0,2.0);

  // draw the true helix
  TCanvas* hcan = new TCanvas("hcan","Traj",1000,1000);
  TPolyLine3D* ttrue = new TPolyLine3D(200);
  TPolyLine3D* tpx = new TPolyLine3D(200);
  TPolyLine3D* tpl = new TPolyLine3D(200);
  TPolyLine3D* tnom = new TPolyLine3D(200);
  Vec4 tpos, xpos, lpos, npos;
  Vec3 tvel, xvel, lvel;
  Mom4 tmom, xmom, lmom;
  float tspeed, xspeed, lspeed;
  double tstep = tptraj.range().range()/200.0;
  for(int istep=0;istep<201;++istep){
  // compute the position from the time
    tpos.SetE(tptraj.range().low() + tstep*istep);
    tptraj.position(tpos);
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
    dxpos->Fill( (xpos-tpos).Vect().R());
    dlpos->Fill( (lpos-tpos).Vect().R());
    tptraj.velocity(tpos.T(),tvel);
    xptraj.velocity(xpos.T(),xvel);
    lptraj.velocity(lpos.T(),lvel);
    dxvel->Fill( (xvel-tvel).R());
    dlvel->Fill( (lvel-tvel).R());
    tspeed = tptraj.speed(tpos.T());
    xspeed = xptraj.speed(tpos.T());
    lspeed = lptraj.speed(tpos.T());
    dxspeed->Fill(xspeed-tspeed);
    dlspeed->Fill(lspeed-tspeed);
    tptraj.momentum(tpos.T(),tmom);
    xptraj.momentum(xpos.T(),xmom);
    lptraj.momentum(lpos.T(),lmom);
    dxmom->Fill( (xmom.Vect()-tmom.Vect()).R());
    dlmom->Fill( (lmom.Vect()-tmom.Vect()).R());
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

  TCanvas* dxcan = new TCanvas("dxcan","dxcan",1200,1200);
  dxcan->Divide(2,2);
  dxcan->cd(1);
  dxpos->Draw();
  dxcan->cd(2);
  dxvel->Draw();
  dxcan->cd(3);
  dxspeed->Draw();
  dxcan->cd(4);
  dxmom->Draw();
  dxcan->Draw();
  dxcan->Write();

  TCanvas* dlcan = new TCanvas("dlcan","dlcan",1200,1200);
  dlcan->Divide(2,2);
  dlcan->cd(1);
  dlpos->Draw();
  dlcan->cd(2);
  dlvel->Draw();
  dlcan->cd(3);
  dlspeed->Draw();
  dlcan->cd(4);
  dlmom->Draw();
  dlcan->Draw();
  dlcan->Write();
 
  return 0;
}


