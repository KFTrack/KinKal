// 
// test derivatives of Traj class
//
#include "KinKal/TPoca.hh"

#include <iostream>
#include <stdio.h>
#include <iostream>
#include <getopt.h>
#include <typeinfo>

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
  printf("Usage: KTrajDerivs  --momentum f --costheta f --azimuth f --particle i --charge i --zorigin f --torigin --dmin f --dmax f --ttest f --By f \n");
}

template <class KTRAJ>
int test(int argc, char **argv) {
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
    {NULL, 0,0,0}


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
      case 'y' : By = atof(optarg);
		 break;
      default: print_usage(); 
	       exit(EXIT_FAILURE);
    }
  }
  // construct original helix from parameters
  Vec3 bnom(0.0,By,1.0);
  Vec4 origin(0.0,0.0,oz,ot);
  float sint = sqrt(1.0-cost*cost);
  // reference helix
  pmass = masses[imass];
  Mom4 momv(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost,pmass);
  KTRAJ refhel(origin,momv,icharge,bnom);
  cout << "Reference " << refhel << endl;
  Vec4 refpos4;
  refpos4.SetE(ttest);
  refhel.position(refpos4);
  cout << "origin position " << origin << " test position " << refpos4 << endl;
  Mom4 refmom;
  refhel.momentum(ttest,refmom);
  int ndel(50);
  // graphs to compare parameter change
  std::vector<TGraph*> pgraphs[3];
  // graphs to compare momentum change
  TGraph* momgraph[3];
  // gaps
  TGraph* gapgraph[3][3];
  // canvases
  TCanvas* dhcan[3];
  TCanvas* dmomcan[3];
 std::string tfname = KTRAJ::trajName() + "Derivs.root";
  TFile lhderiv(tfname.c_str(),"RECREATE");
  // loop over derivative directions
  double del = (dmax-dmin)/(ndel-1);
  for(int idir=0;idir<3;++idir){
    KInter::MDir tdir =static_cast<KInter::MDir>(idir);
//    cout << "testing direction " << KInter::directionName(tdir) << endl;
    // parameter change
    pgraphs[idir] = std::vector<TGraph*>(KTRAJ::NParams(),0); 
    for(size_t ipar = 0; ipar < KTRAJ::NParams(); ipar++){
      pgraphs[idir][ipar] = new TGraph(ndel);
      string title = KTRAJ::paramName(typename KTRAJ::ParamIndex(ipar));
      title += ";exact;1st derivative";
      pgraphs[idir][ipar]->SetTitle(title.c_str());
    }
    momgraph[idir] = new TGraph(ndel);
    momgraph[idir]->SetTitle("Momentum Direction;exact;1st derivative");
    for(int jdir=0;jdir < 3;jdir++){
      KInter::MDir tjdir =static_cast<KInter::MDir>(jdir);
      gapgraph[idir][jdir] = new TGraph(ndel);
      string title = "Gap in " + KInter::directionName(tjdir) + ";change;gap value (mm)";
      gapgraph[idir][jdir]->SetTitle(title.c_str());
    }
    // scan range of change
    for(int id=0;id<ndel;++id){
      double delta = dmin + del*id; 
//      cout << "Delta = " << delta << endl;
      // compute 1st order change in parameters
      typename KTRAJ::DVEC pder;
      Vec3 dmomdir;
      refhel.momDeriv(tdir,ttest,pder,dmomdir);
      //  compute exact altered params
      Vec3 newmom = refmom.Vect() + delta*dmomdir*mom;
      Mom4 momv(newmom.X(),newmom.Y(),newmom.Z(),pmass);
      KTRAJ xhel(refpos4,momv,icharge,bnom);
//      cout << "derivative vector" << pder << endl;
      auto dvec = refhel.params().parameters() + delta*pder;
      typename KTRAJ::PDATA pdata(dvec,refhel.params().covariance());
      KTRAJ dhel(pdata,refhel.mass(),refhel.charge(),bnom);
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
      // project along 3 directions
      for(int jdir=0;jdir < 3;jdir++){
	KInter::MDir tjdir =static_cast<KInter::MDir>(jdir);
	Vec3 jmomdir;
	refhel.momDeriv(tjdir,ttest,pder,jmomdir);
	gapgraph[idir][jdir]->SetPoint(id,delta,gap.Vect().Dot(jmomdir));
      }
      // parameter diff
      for(size_t ipar = 0; ipar < KTRAJ::NParams(); ipar++){
	pgraphs[idir][ipar]->SetPoint(id,xhel.paramVal(ipar)-refhel.paramVal(ipar),dhel.paramVal(ipar)-refhel.paramVal(ipar));
      }
      // compare momenta after change
      //
      Vec3 dxmom = momv.Vect() - refmom.Vect();
      Vec3 ddmom = dmom.Vect() - refmom.Vect();
      momgraph[idir]->SetPoint(id,dxmom.Dot(dmomdir),ddmom.Dot(dmomdir));
    }
    char title[80];
    char name[80];
    snprintf(name,80,"dh%s",KInter::directionName(tdir).c_str());
    snprintf(title,80,"Helix Change %s",KInter::directionName(tdir).c_str());
    dhcan[idir] = new TCanvas(name,title,1200,800);
    dhcan[idir]->Divide(3,2);
    for(size_t ipar = 0; ipar < KTRAJ::NParams(); ipar++){
      dhcan[idir]->cd(ipar+1);
      pgraphs[idir][ipar]->Draw("AC*");
    }
    dhcan[idir]->Draw();
    dhcan[idir]->Write();

    snprintf(name,80,"dm%s",KInter::directionName(tdir).c_str());
    snprintf(title,80,"Mom Change %s",KInter::directionName(tdir).c_str());
    dmomcan[idir] = new TCanvas(name,title,800,800);
    dmomcan[idir]->Divide(2,2);
    dmomcan[idir]->cd(1);
    momgraph[idir]->Draw("AC*");
    for(int jdir=0;jdir < 3;jdir++){
      dmomcan[idir]->cd(2+jdir);
      gapgraph[idir][jdir]->Draw("AC*");
    }
    dmomcan[idir]->Draw();
    dmomcan[idir]->Write();
  }

  lhderiv.Write();
  lhderiv.Close();
  return 0;
}

