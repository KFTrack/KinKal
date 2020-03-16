// 
// test derivatives of the Loop Helix TTraj class
//
#include "KinKal/LHelix.hh"
#include "KinKal/Context.hh"

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
  printf("Usage: LHelixDerivs  --momentum f --costheta f --azimuth f --particle i --charge i --zorigin f --torigin --dmin f --dmax f --ttest f\n");
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
    {"ttest",     required_argument, 0, 't'  }


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
      default: print_usage(); 
	       exit(EXIT_FAILURE);
    }
  }
  // construct original helix from parameters
  UniformBField BF(1.0);
  Context context(BF);
  Vec4 origin(0.0,0.0,oz,ot);
  float sint = sqrt(1.0-cost*cost);
  // reference helix
  pmass = masses[imass];
  Mom4 momv(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost,pmass);
  LHelix refhel(origin,momv,icharge,context);
  cout << "Reference " << refhel << endl;
  Vec4 refpos4;
  refpos4.SetE(ttest);
  refhel.position(refpos4);
  cout << "origin position " << origin << " test position " << refpos4 << endl;
  Mom4 refmom;
  refhel.momentum(ttest,refmom);
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
  // position change
  TGraph* posgraph[6];
  // gaps
  TGraph* gapgraph[3];
  // canvases
  TCanvas* dhcan[3];
  TCanvas* dmomcan[3];
  TFile lhderiv("LHelixDerivs.root","RECREATE");
  // loop over derivative directions
  double del = (dmax-dmin)/(ndel-1);
  for(int idir=0;idir<3;++idir){
    KInter::MDir tdir =static_cast<KInter::MDir>(idir);
    Vec3 dmomdir;
    refhel.dirVector(tdir,ttest,dmomdir);
//    cout << "testing direction " << KInter::directionName(tdir) << endl;
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
    for(int id=0;id<ndel;++id){
      double delta = dmin + del*id; 
//      cout << "Delta = " << delta << endl;
      //  compute exact altered params
      Vec3 newmom = refmom.Vect() + delta*dmomdir*mom;
      Mom4 momv(newmom.X(),newmom.Y(),newmom.Z(),pmass);
      LHelix xhel(refpos4,momv,icharge,context);
      // now, compute 1st order change in parameters
      LHelix::PDER pder;
      refhel.momDeriv(tdir,ttest,pder);
//      cout << "derivative vector" << pder << endl;
      auto dvec = refhel.params().parameters() + delta*pder;
      LHelix dhel(dvec,refhel.params().covariance(),refhel.mass(),refhel.charge(),context);
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
      Vec4 gap = dpos - refpos4;
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
      refhel.dirVector(KInter::momdir,ttest,changedir);
      mom0graph[idir]->SetPoint(id,dxmom.Dot(changedir),ddmom.Dot(changedir));
      refhel.dirVector(KInter::theta1,ttest,changedir);
      mom1graph[idir]->SetPoint(id,dxmom.Dot(changedir),ddmom.Dot(changedir));
      refhel.dirVector(KInter::theta2,ttest,changedir);
      mom2graph[idir]->SetPoint(id,dxmom.Dot(changedir),ddmom.Dot(changedir));
    }
    char title[80];
    char name[80];
    snprintf(name,80,"dh%s",KInter::directionName(tdir).c_str());
    snprintf(title,80,"Helix Change %s",KInter::directionName(tdir).c_str());
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
    dhcan[idir]->Write();

    snprintf(name,80,"dm%s",KInter::directionName(tdir).c_str());
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
    dmomcan[idir]->Write();
  }

// now spatial derivatives: these are used to constrain continuity between traj pieces
  LHelix::PDER pder;
  refhel.posDeriv(ttest,pder);
  Vec3 refpos, refdir;
  refhel.position(ttest,refpos);
  refhel.direction(ttest,refdir);
  double refpdot = refpos.Dot(refdir);

  TCanvas* dpdotcan = new TCanvas("dpdot","P dot D derivatives",1200,800);
  dpdotcan->Divide(3,2);
  for(size_t ipar=0;ipar< LHelix::NParams();ipar++){
    LHelix::ParamIndex ip = static_cast<LHelix::ParamIndex>(ipar);
    posgraph[ipar] = new TGraph(ndel);
    string title = LHelix::paramName(ip) + string(" #Delta #vec{P} #bullet #vec{D};exact;1st derivative");
    posgraph[ipar]->SetTitle(title.c_str());
    LHelix xhel(refhel);
    for(int id=0;id<ndel;++id){
      double delta = dmin + del*id;
      xhel.params().parameters()[ipar] = refhel.params().parameters()[ipar]*(1.0 + delta);
      Vec3 pos;
      xhel.position(ttest,pos);
      double dpdot = pos.Dot(refdir) - refpdot;
      // linear approximation
      double dirdpdot = pder[ipar]*delta*refhel.params().parameters()[ipar];
      posgraph[ipar]->SetPoint(id,dpdot,dirdpdot);
    }
    dpdotcan->cd(ipar+1);
    posgraph[ipar]->Draw("AC*");
  }
  dpdotcan->Write();

  lhderiv.Write();
  lhderiv.Close();
  return 0;
}

