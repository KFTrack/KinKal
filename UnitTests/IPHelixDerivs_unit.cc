// 
// test derivatives of the Loop Helix TTraj class
//
#include "KinKal/IPHelix.hh"

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
#include "TDirectory.h"
#include "TProfile.h"
#include "TProfile2D.h"

using namespace KinKal;
using namespace std;

void print_usage() {
  printf("Usage: IPHelixDerivs  --momentum f --costheta f --azimuth f --particle i --charge i --zorigin f --torigin --dmin f --dmax f --ttest f\n");
}

int main(int argc, char **argv) {
  gROOT->SetBatch(kTRUE);
  // save canvases
  int opt;
  double mom(105.0), cost(0.7), phi(0.5);
  double masses[5]={0.511,105.66,139.57, 493.68, 938.0};
  int imass(0), icharge(-1);
  double pmass, oz(100.0), ot(0.0), ttest(5.0);
  double dmin(-5e-2), dmax(5e-2);
  double By(0.0);

  static struct option long_options[] = {
    {"momentum",     required_argument, 0, 'm' },
    {"costheta",     required_argument, 0, 'c'  },
    {"azimuth",     required_argument, 0, 'a'  },
    {"particle",     required_argument, 0, 'p'  },
    {"charge",     required_argument, 0, 'q'  },
    {"zorigin",     required_argument, 0, 'z'  },
    {"torigin",     required_argument, 0, 'o'  },
    {"dmin",     required_argument, 0, 's'  },
    {"dmax",     required_argument, 0, 'e'  },
    {"ttest",     required_argument, 0, 't'  },
    {"By",     required_argument, 0, 'y'  },


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
      case 'o' : ot = atof(optarg);
		 break;
      case 's' : dmin = atof(optarg);
		 break;
      case 'e' : dmax = atof(optarg);
		 break;
      case 't' : ttest = atof(optarg);
    break;
      case 'y': By = atof(optarg);
    break;
      default: print_usage();
	       exit(EXIT_FAILURE);
    }
  }
  // construct original helix from parameters
  Vec3 bnom(0.0, By, 1.0);
  Vec4 origin(0.0,0.0,oz,ot);
  float sint = sqrt(1.0-cost*cost);
  // reference helix
  pmass = masses[imass];
  Mom4 momv(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost,pmass);
  IPHelix refhel(origin,momv,icharge,bnom);
  cout << "Reference " << refhel << endl;
  Vec4 refpos;
  refpos.SetE(ttest);
  refhel.position(refpos);
  cout << "origin position " << origin << " test position " << refpos << endl;
  Mom4 refmom;
  refhel.momentum(ttest,refmom);
  int ndel(50);
  // graphs to compare parameter change
  TGraph* d0graph[3];
  TGraph* omegagraph[3];
  TGraph* t0graph[3];
  TGraph* phi0graph[3];
  TGraph* tanDipgraph[3];
  TGraph* z0graph[3];
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
    KInter::MDir tdir = static_cast<KInter::MDir>(idir);
    Vec3 dmomdir;
    refhel.dirVector(tdir,ttest,dmomdir);
//    cout << "testing direction " << KInter::directionName(tdir) << endl;
    // parameter change

    d0graph[idir] = new TGraph(ndel);
    d0graph[idir]->SetTitle("d_{0};exact;1st derivative");
    omegagraph[idir] = new TGraph(ndel);
    omegagraph[idir]->SetTitle("#omega;exact;1st derivative");
    t0graph[idir] = new TGraph(ndel);
    t0graph[idir]->SetTitle("t_{0};exact;1st derivative");
    phi0graph[idir] = new TGraph(ndel);
    phi0graph[idir]->SetTitle("#phi_{0};exact;1st derivative");
    tanDipgraph[idir] = new TGraph(ndel);
    tanDipgraph[idir]->SetTitle("tan#lambda;exact;1st derivative");
    z0graph[idir] = new TGraph(ndel);
    z0graph[idir]->SetTitle("z_{0};exact;1st derivative");
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
    for(int id=0;id<ndel;++id){
      double delta = dmin + del*id; 
//      cout << "Delta = " << delta << endl;
      //  compute exact altered params
      Vec3 newmom = refmom.Vect() + delta*dmomdir*mom;
      Mom4 momv(newmom.X(),newmom.Y(),newmom.Z(),pmass);
      IPHelix xhel(refpos,momv,icharge,bnom);
      // now, compute 1st order change in parameters
      IPHelix::PDER pder;
      refhel.momDeriv(tdir,ttest,pder);
//      cout << "derivative vector" << pder << endl;
      auto dvec = refhel.params().parameters();
      for (size_t ipar = 0; ipar < 6; ipar++)
        dvec[ipar] += delta * pder[ipar];
      //dvec += pder;
      IPHelix dhel(dvec, refhel.params().covariance(), refhel.mass(), refhel.charge(), bnom);
      // test
      Vec4 xpos, dpos;
      xpos.SetE(ttest);
      dpos.SetE(ttest);
      xhel.position(xpos);
      dhel.position(dpos);
//      cout << " exa pos " << xpos << endl
//      << " del pos " << dpos << endl;
      Mom4 dmom;
      dhel.momentum(ttest,dmom);
//      cout << "Exact change" << xhel << endl;
//      cout << "Derivative  " << dhel << endl;
      Vec4 gap = dpos - refpos;
      gapgraph[idir]->SetPoint(id,delta,sqrt(gap.Vect().Mag2()));
      // parameter diff
      d0graph[idir]->SetPoint(id,xhel.d0()-refhel.d0(),dhel.d0()-refhel.d0());
      omegagraph[idir]->SetPoint(id,xhel.omega()-refhel.omega(),dhel.omega()-refhel.omega());
      t0graph[idir]->SetPoint(id,xhel.t0()-refhel.t0(),dhel.t0()-refhel.t0());
      phi0graph[idir]->SetPoint(id,xhel.phi0()-refhel.phi0(),dhel.phi0()-refhel.phi0());
      tanDipgraph[idir]->SetPoint(id,xhel.tanDip()-refhel.tanDip(),dhel.tanDip()-refhel.tanDip());
      z0graph[idir]->SetPoint(id,xhel.z0()-refhel.z0(),dhel.z0()-refhel.z0());

      // compare momenta after change
      //
      Vec3 dxmom = momv.Vect() - refmom.Vect();
      Vec3 ddmom = dmom.Vect() - refmom.Vect();
      Vec3 changedir;
      refhel.dirVector(KInter::momdir, ttest, changedir);
      mom0graph[idir]->SetPoint(id, dxmom.Dot(changedir), ddmom.Dot(changedir));
      refhel.dirVector(KInter::theta1, ttest, changedir);
      mom1graph[idir]->SetPoint(id, dxmom.Dot(changedir), ddmom.Dot(changedir));
      refhel.dirVector(KInter::theta2, ttest, changedir);
      mom2graph[idir]->SetPoint(id, dxmom.Dot(changedir), ddmom.Dot(changedir));
    }

    // draw comparisons
    char title[80];
    char name[80];
    snprintf(name,80,"dhcan%s",KInter::directionName(tdir).c_str());
    snprintf(title,80,"Helix Change %s",KInter::directionName(tdir).c_str());
    dhcan[idir] = new TCanvas(name,title,1200,800);
    dhcan[idir]->Divide(3,2);
    dhcan[idir]->cd(1);
    d0graph[idir]->Draw("AC*");
    dhcan[idir]->cd(2);
    omegagraph[idir]->Draw("AC*");
    dhcan[idir]->cd(3);
    t0graph[idir]->Draw("AC*");
    dhcan[idir]->cd(4);
    phi0graph[idir]->Draw("AC*");
    dhcan[idir]->cd(5);
    tanDipgraph[idir]->Draw("AC*");
    dhcan[idir]->cd(6);
    z0graph[idir]->Draw("AC*");
    dhcan[idir]->Draw();

    snprintf(name,80,"dmcan_%s",KInter::directionName(tdir).c_str());
    snprintf(title,80,"Mom Change %s",KInter::directionName(tdir).c_str());
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
    snprintf(fname,100,"IPHelixDerivs_dh_%s.root",KInter::directionName(tdir).c_str());
    dhcan[idir]->SaveAs(fname);
    snprintf(fname,100,"IPHelixDerivs_dmom_%s.root",KInter::directionName(tdir).c_str());
    dmomcan[idir]->SaveAs(fname);

  }

 return 0;
}

