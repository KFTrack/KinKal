// 
// test derivatives of the Loop Helix TTraj class
//
#include "BTrk/KinKal/LHelix.hh"
#include "BTrk/KinKal/Context.hh"

#include <iostream>
#include <stdio.h>
#include <iostream>
#include <getopt.h>

#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
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
  printf("Usage: LHelixDerivs  --momentum f --costheta f --azimuth f --particle i --charge i --zorigin f --torigin --dmin ff--dmax f \n");
}

int main(int argc, char **argv) {
  gROOT->SetBatch(kTRUE);
  // save canvases
  int opt;
  double mom(105.0), cost(0.7), phi(0.5);
  double masses[5]={0.511,105.66,139.57, 493.68, 938.0};
  int imass(0), icharge(-1);
  double pmass, oz(100.0), ot(0.0);
  double dmin(-5e-2), dmax(5e-2);

  static struct option long_options[] = {
    {"momentum",     required_argument, 0, 'm' },
    {"costheta",     required_argument, 0, 'c'  },
    {"azimuth",     required_argument, 0, 'a'  },
    {"particle",     required_argument, 0, 'p'  },
    {"charge",     required_argument, 0, 'q'  },
    {"zorigin",     required_argument, 0, 'z'  },
    {"torigin",     required_argument, 0, 't'  },
    {"tmin",     required_argument, 0, 's'  },
    {"tmax",     required_argument, 0, 'e'  }


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
      case 'z' : oz = atof(optarg);
		 break;
      case 't' : ot = atof(optarg);
		 break;
      case 's' : dmin = atof(optarg);
		 break;
      case 'e' : dmax = atof(optarg);
		 break;
      default: print_usage(); 
	       exit(EXIT_FAILURE);
    }
  }
  // construct original helix from parameters
  Context context;
  context.Bz_ = 1.0; // 1 Tesla
  Vec4 origin(0.0,0.0,oz,ot);
  float sint = sqrt(1.0-cost*cost);
  // reference helix
  pmass = masses[imass];
  Mom4 momv(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost,pmass);
  LHelix refhel(origin,momv,icharge,context);
  cout << "Reference " << refhel << endl;
  Vec4 refpos;
  refpos.SetE(ot);
  refhel.position(refpos);
  cout << "origin " << origin << " reference " << refpos << endl;
  Mom4 refmom;
  refhel.momentum(ot,refmom);
  int ndel(50);
  // graphs to compare parameter change
  TGraph* radgraph[3];
  TGraph* lambdagraph[3];
  TGraph* t0graph[3];
  TGraph* phi0graph[3];
  TGraph* cxgraph[3];
  TGraph* cygraph[3];
  // graphs to compare momentum change
  TGraph* mom0graph[3];
  TGraph* mom1graph[3];
  TGraph* mom2graph[3];
  // gaps
  TGraph* gapgraph[3];
  // canvases
  TCanvas* dhcan[3];
  TCanvas* dmomcan[3];

  // loop over derivative directions
  for(int idir=0;idir<3;++idir){
    KTraj::trajdir tdir =static_cast<KTraj::trajdir>(idir);
//    cout << "testing direction " << KTraj::directionName(tdir) << endl;
    // parameter change

    radgraph[idir] = new TGraph(ndel);
    radgraph[idir]->SetTitle("Radius;exact;1st derivative");
    lambdagraph[idir] = new TGraph(ndel);
    lambdagraph[idir]->SetTitle("Lambda;exact;1st derivative");
    t0graph[idir] = new TGraph(ndel);
    t0graph[idir]->SetTitle("t_{0};exact;1st derivative");
    phi0graph[idir] = new TGraph(ndel);
    phi0graph[idir]->SetTitle("phi_{0};exact;1st derivative");
    cxgraph[idir] = new TGraph(ndel);
    cxgraph[idir]->SetTitle("C_{x};exact;1st derivative");
    cygraph[idir] = new TGraph(ndel);
    cygraph[idir]->SetTitle("C_{y};exact;1st derivative");
    mom0graph[idir] = new TGraph(ndel);
    mom0graph[idir]->SetTitle("Momentum Direction;exact;1st derivative");
    mom1graph[idir] = new TGraph(ndel);
    mom1graph[idir]->SetTitle("Theta Direction;exact;1st derivative");
    mom2graph[idir] = new TGraph(ndel);
    mom2graph[idir]->SetTitle("Phi Direction;exact;1st derivative");
    gapgraph[idir] = new TGraph(ndel);
    gapgraph[idir]->SetTitle("Gap;change;gap value (mm)");

    // scan range of change
    double del = (dmax-dmin)/(ndel-1);
    Mom4 momv;
    for(int id=0;id<ndel;++id){
      double delta = dmin + del*id; 
//      cout << "Delta = " << delta << endl;
 
      //  compute exact altered params
      if(tdir == KTraj::theta1){
	double newcost = cos(acos(cost) + delta);
	double newsint = sqrt(1.0-newcost*newcost);
	momv = Mom4(mom*newsint*cos(phi),mom*newsint*sin(phi),mom*newcost,pmass);
      } else if(tdir == KTraj::theta2){
	double dphi = delta/sint;
	double newphi = phi + dphi;
	momv = Mom4(mom*sint*cos(newphi),mom*sint*sin(newphi),mom*cost,pmass);
      } else if(tdir == KTraj::momdir) {
	double newmom = (1.0+delta)*mom;
	momv = Mom4(newmom*sint*cos(phi),newmom*sint*sin(phi),newmom*cost,pmass);
      }
      LHelix xhel(origin,momv,icharge,context);
      // now, compute 1st order change in parameters
      LHelix::PDer pder;
      refhel.momDeriv(tdir,ot,pder);
//      cout << "derivative vector" << pder << endl;
      LHelix::TDATA::DVec dvec = refhel.params().vec();
      for(size_t ipar=0;ipar<6;ipar++)
	dvec[ipar] += delta*pder[ipar][0];
      //dvec += pder;
      LHelix dhel(dvec,refhel.params().mat(),refhel.mass(),refhel.charge(),context);
      Mom4 dmom;
      dhel.momentum(ot,dmom);
//      cout << "Exact change" << xhel << endl;
//      cout << "Derivative  " << dhel << endl;
      Vec4 dpos;
      dpos.SetE(ot);
      dhel.position(dpos);
      Vec4 gap = dpos - refpos;
      gapgraph[idir]->SetPoint(id,delta,sqrt(gap.Vect().Mag2()));
      // parameter diff
      radgraph[idir]->SetPoint(id,xhel.rad()-refhel.rad(),dhel.rad()-refhel.rad());
      lambdagraph[idir]->SetPoint(id,xhel.lam()-refhel.lam(),dhel.lam()-refhel.lam());
      t0graph[idir]->SetPoint(id,xhel.t0()-refhel.t0(),dhel.t0()-refhel.t0());
      phi0graph[idir]->SetPoint(id,xhel.phi0()-refhel.phi0(),dhel.phi0()-refhel.phi0());
      cxgraph[idir]->SetPoint(id,xhel.cx()-refhel.cx(),dhel.cx()-refhel.cx());
      cygraph[idir]->SetPoint(id,xhel.cy()-refhel.cy(),dhel.cy()-refhel.cy());

      // compare momenta after change
      //
      Vec3 dxmom = momv.Vect() - refmom.Vect();
      Vec3 ddmom = dmom.Vect() - refmom.Vect();
      Vec3 changedir;
      refhel.dirVector(KTraj::momdir,ot,changedir);
      mom0graph[idir]->SetPoint(id,dxmom.Dot(changedir),ddmom.Dot(changedir));
      refhel.dirVector(KTraj::theta1,ot,changedir);
      mom1graph[idir]->SetPoint(id,dxmom.Dot(changedir),ddmom.Dot(changedir));
      refhel.dirVector(KTraj::theta2,ot,changedir);
      mom2graph[idir]->SetPoint(id,dxmom.Dot(changedir),ddmom.Dot(changedir));
    }
    // draw comparisons
    char title[80];
    char name[80];
    snprintf(name,80,"dhcan%s",KTraj::directionName(tdir).c_str());
    snprintf(title,80,"Helix Change %s",KTraj::directionName(tdir).c_str());
    dhcan[idir] = new TCanvas(name,title,1200,800);
    dhcan[idir]->Divide(3,2);
    dhcan[idir]->cd(1);
    radgraph[idir]->Draw("AC*");
    dhcan[idir]->cd(2);
    lambdagraph[idir]->Draw("AC*");
    dhcan[idir]->cd(3);
    t0graph[idir]->Draw("AC*");
    dhcan[idir]->cd(4);
    phi0graph[idir]->Draw("AC*");
    dhcan[idir]->cd(5);
    cxgraph[idir]->Draw("AC*");
    dhcan[idir]->cd(6);
    cygraph[idir]->Draw("AC*");
    dhcan[idir]->Draw();

    snprintf(name,80,"dmcan_%s",KTraj::directionName(tdir).c_str());
    snprintf(title,80,"Mom Change %s",KTraj::directionName(tdir).c_str());
    dmomcan[idir] = new TCanvas(name,title,800,800);
    dmomcan[idir]->Divide(2,2);
    dmomcan[idir]->cd(1);
    mom0graph[idir]->Draw("AC*");
    dmomcan[idir]->cd(2);
    mom1graph[idir]->Draw("AC*");
    dmomcan[idir]->cd(3);
    mom2graph[idir]->Draw("AC*");
    dmomcan[idir]->cd(4);
    gapgraph[idir]->Draw("AC*");
    dmomcan[idir]->Draw();
  
    char fname[100];
    snprintf(fname,100,"LHelixDerivs_dh_%s.root",KTraj::directionName(tdir).c_str());
    dhcan[idir]->SaveAs(fname);
    snprintf(fname,100,"LHelixDerivs_dmom_%s.root",KTraj::directionName(tdir).c_str());
    dmomcan[idir]->SaveAs(fname);

  }

 return 0;
}

