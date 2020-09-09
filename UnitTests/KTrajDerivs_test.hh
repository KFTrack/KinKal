// 
// test derivatives of Traj class
//
#include "KinKal/ClosestApproach.hh"

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
#include "TF1.h"
#include "TRandom3.h"
#include "TH2F.h"
#include "TDirectory.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TFitResult.h"

using namespace KinKal;
using namespace std;

void print_usage() {
  printf("Usage: KTrajDerivs  --momentum f --costheta f --azimuth f --particle i --charge i --zorigin f --torigin --delta f --ttest f --By f \n");
}

template <class KTRAJ>
int test(int argc, char **argv) {
  gROOT->SetBatch(kTRUE);
  // save canvases
  int status(0);
  int opt;
  double mom(105.0), cost(0.7), phi(0.5);
  double masses[5]={0.511,105.66,139.57, 493.68, 938.0};
  int imass(0), icharge(-1);
  double pmass, oz(100.0), ot(0.0), ttest(5.0);
  double delta(1e-2);
  double By(0.0);

  static struct option long_options[] = {
    {"momentum",     required_argument, 0, 'm' },
    {"costheta",     required_argument, 0, 'c'  },
    {"azimuth",     required_argument, 0, 'a'  },
    {"particle",     required_argument, 0, 'p'  },
    {"charge",     required_argument, 0, 'q'  },
    {"zorigin",     required_argument, 0, 'z'  },
    {"torigin",     required_argument, 0, 'o'  },
    {"delta",     required_argument, 0, 'd'  },
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
      case 'd' : delta = atof(optarg);
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
  VEC3 bnom(0.0,By,1.0);
  VEC4 origin(0.0,0.0,oz,ot);
  double sint = sqrt(1.0-cost*cost);
  // reference helix
  pmass = masses[imass];
  MOM4 momv(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost,pmass);
  KTRAJ refhel(origin,momv,icharge,bnom);
  //cout << "Reference " << refhel << endl;
  VEC4 refpos4;
  refpos4.SetE(ttest);
  refhel.position(refpos4);
  //  cout << "origin position " << origin << " test position " << refpos4 << endl;
  MOM4 refmom = refhel.momentum(ttest);
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
  double del = 2*delta/(ndel-1);
  double dmin = -delta;
  for(int idir=0;idir<3;++idir){
    MomBasis::Direction tdir =static_cast<MomBasis::Direction>(idir);
//    cout << "testing direction " << MomBasis::directionName(tdir) << endl;
    // parameter change
    pgraphs[idir] = std::vector<TGraph*>(NParams(),0); 
    for(size_t ipar = 0; ipar < NParams(); ipar++){
      pgraphs[idir][ipar] = new TGraph(ndel);
      string title = KTRAJ::paramName(typename KTRAJ::ParamIndex(ipar));
      title += ";exact;1st derivative";
      pgraphs[idir][ipar]->SetTitle(title.c_str());
    }
    momgraph[idir] = new TGraph(ndel);
    momgraph[idir]->SetTitle("Momentum Direction;exact;1st derivative");
    for(int jdir=0;jdir < 3;jdir++){
      MomBasis::Direction tjdir =static_cast<MomBasis::Direction>(jdir);
      gapgraph[idir][jdir] = new TGraph(ndel);
      string title = "Gap in " + MomBasis::directionName(tjdir) + ";Fractional change;Gap value (mm)";
      gapgraph[idir][jdir]->SetTitle(title.c_str());
    }
    // scan range of change
    DVEC pder = refhel.momDeriv(ttest,tdir);
    for(int id=0;id<ndel;++id){
      double dval = dmin + del*id;
//      cout << "Delta = " << dval << endl;
      // compute 1st order change in parameters
      VEC3 dmomdir = refhel.direction(ttest,tdir);
      //  compute exact altered params
      VEC3 newmom = refmom.Vect() + dval*dmomdir*mom;
      MOM4 momv(newmom.X(),newmom.Y(),newmom.Z(),pmass);
      KTRAJ xhel(refpos4,momv,icharge,bnom);
//      cout << "derivative vector" << pder << endl;
      DVEC dvec = refhel.params().parameters() + dval*pder;
      Parameters pdata(dvec,refhel.params().covariance());
      KTRAJ dhel(pdata,refhel);
      // test
      VEC4 xpos, dpos;
      xpos.SetE(ttest);
      dpos.SetE(ttest);
      xhel.position(xpos);
      dhel.position(dpos);
//      cout << " exa pos " << xpos << endl
//      << " del pos " << dpos << endl;
      MOM4 dmom = dhel.momentum(ttest);
//      cout << "Exact change" << xhel << endl;
//      cout << "Derivative  " << dhel << endl;
      VEC4 gap = dpos - refpos4;
      // project along 3 directions
      for(int jdir=0;jdir < 3;jdir++){
	MomBasis::Direction tjdir =static_cast<MomBasis::Direction>(jdir);
	VEC3 jmomdir = refhel.direction(ttest,tjdir);
	gapgraph[idir][jdir]->SetPoint(id,dval,gap.Vect().Dot(jmomdir));
      }
      // parameter diff
      for(size_t ipar = 0; ipar < NParams(); ipar++){
	pgraphs[idir][ipar]->SetPoint(id,xhel.paramVal(ipar)-refhel.paramVal(ipar),dhel.paramVal(ipar)-refhel.paramVal(ipar));
      }
      // compare momenta after change
      //
      VEC3 dxmom = momv.Vect() - refmom.Vect();
      VEC3 ddmom = dmom.Vect() - refmom.Vect();
      momgraph[idir]->SetPoint(id,dxmom.Dot(dmomdir),ddmom.Dot(dmomdir));
    }
    char gtitle[80];
    char gname[80];
    snprintf(gname,80,"dh%s",MomBasis::directionName(tdir).c_str());
    snprintf(gtitle,80,"KTraj Change %s",MomBasis::directionName(tdir).c_str());
    dhcan[idir] = new TCanvas(gname,gtitle,1200,800);
    dhcan[idir]->Divide(3,2);
    TF1* pline = new TF1("pline","[0]+[1]*x");
    for(size_t ipar = 0; ipar < NParams(); ipar++){
      dhcan[idir]->cd(ipar+1);
      // if this is non-trivial, fit
      if(fabs(pder[ipar])>1e-9){
	pline->SetParameters(0.0,1.0);
	TFitResultPtr pfitr = pgraphs[idir][ipar]->Fit(pline,"SQ","AC*");
	pgraphs[idir][ipar]->Draw("AC*");
	if(fabs(pfitr->Parameter(0))> 10*delta || fabs(pfitr->Parameter(1)-1.0) > delta){
	  cout << "Parameter " 
	    << KTRAJ::paramName(typename KTRAJ::ParamIndex(ipar))
	    << " in direction " << MomBasis::directionName(tdir)
	    << " Out of tolerance : Offset " << pfitr->Parameter(0) << " Slope " << pfitr->Parameter(1) << endl;
	  status = 1;
	}
      }
    }
    dhcan[idir]->Draw();
    dhcan[idir]->Write();

    snprintf(gname,80,"dm%s",MomBasis::directionName(tdir).c_str());
    snprintf(gtitle,80,"Mom Change %s",MomBasis::directionName(tdir).c_str());
    dmomcan[idir] = new TCanvas(gname,gtitle,800,800);
    dmomcan[idir]->Divide(2,2);
    dmomcan[idir]->cd(1);
    pline->SetParameters(0.0,1.0);
    TFitResultPtr pfitr = momgraph[idir]->Fit(pline,"SQ","AC*");
    momgraph[idir]->Draw("AC*");
    if(fabs(pfitr->Parameter(0))> 10*delta || fabs(pfitr->Parameter(1)-1.0) > 0.1*delta){
      cout << "Momentum Direction " 
	<< MomBasis::directionName(tdir)
	<< " Out of tolerance : Offset " << pfitr->Parameter(0) << " Slope " << pfitr->Parameter(1) << endl;
      status = 1;
    }
    for(int jdir=0;jdir < 3;jdir++){
      dmomcan[idir]->cd(2+jdir);
      gapgraph[idir][jdir]->Draw("AC*");

    }
    dmomcan[idir]->Draw();
    dmomcan[idir]->Write();
  }

  // test parameter<->phase space translation
  auto dPdX = refhel.dPardX(ttest);  
  auto dPdM = refhel.dPardM(ttest);  
  auto dXdP = refhel.dXdPar(ttest);  
  auto dMdP = refhel.dMdPar(ttest);  
  auto dPdS = refhel.dPardState(ttest);
  auto dSdP = refhel.dStatedPar(ttest);
  auto ptest = dPdS*dSdP;
  for(size_t irow=0;irow<NParams();irow++) {
    for(size_t icol=0;icol<NParams();icol++) {
      double val(0.0);
      if(irow==icol)val = 1.0;
      if(fabs(ptest(irow,icol) - val) > 1e-9){
	cout <<"Error in parameter derivative test" << endl;
	status = 1;
      }
    }
  }
  auto xtest = dXdP*dPdX;
  for(size_t irow=0;irow<3;irow++) {
    for(size_t icol=0;icol<3;icol++) {
      double val(0.0);
      if(irow==icol)val = 1.0;
      if(fabs(xtest(irow,icol) - val) > 1e-9){
	cout <<"Error in position derivative test" << endl;
	status = 1;
      }
    }
  }
  auto mtest = dMdP*dPdM;
  for(size_t irow=0;irow<3;irow++) {
    for(size_t icol=0;icol<3;icol++) {
      double val(0.0);
      if(irow==icol)val = 1.0;
      if(fabs(mtest(irow,icol) - val) > 1e-9){
	cout <<"Error in momentum derivative test" << endl;
	status = 1;
      }
    }
  }

// test changes due to BFieldMap
  TCanvas* dbcan[3]; // 3 directions
  std::vector<TGraph*> bpgraphs[3];
  std::array<VEC3,3> basis = {VEC3(1.0,0.0,0.0), VEC3(0.0,1.0,0.0), VEC3(0.0,0.0,1.0) };
  std::array<std::string,3> anames = {"X", "Y", "Z"};
  // gaps
  TGraph* bgapgraph[3];
  for(int idir=0;idir<3;++idir){
    bpgraphs[idir] = std::vector<TGraph*>(NParams(),0); 
    for(size_t ipar = 0; ipar < NParams(); ipar++){
      bpgraphs[idir][ipar] = new TGraph(ndel);
      string title = KTRAJ::paramName(typename KTRAJ::ParamIndex(ipar));
      title += ";exact;1st derivative";
      bpgraphs[idir][ipar]->SetTitle(title.c_str());
    }
    bgapgraph[idir] = new TGraph(ndel);
    string title = "Gap for #Delta B in " + anames[idir] + ";Fractional change;Gap value (mm)";
    bgapgraph[idir]->SetTitle(title.c_str());
    for(int id=0;id<ndel;++id){
      // construct exact helix for this field and the corresponding exact parameter change
      double dval = dmin + del*id;
      VEC3 bf = bnom + basis[idir]*dval;
      auto state = refhel.measurementState(ttest);
      // exact traj given the full state
      KTRAJ newbfhel(state,ttest,refhel.mass(),refhel.charge(),bf);
      auto newstate = newbfhel.measurementState(ttest);
      for(size_t ipar=0;ipar < ParticleState::dimension(); ipar++){
	if(fabs(state.stateVector().state()[ipar] - newstate.stateVector().state()[ipar])>1.0e-6) cout << "State vector " << ipar << " doesn't match: original "
	  << state.stateVector().state()[ipar] << " rotated " << newstate.stateVector().state()[ipar]  << endl;
      }
      DVEC dpx = newbfhel.params().parameters() - refhel.params().parameters();
      // 1st order change trajectory
      KTRAJ dbtraj(refhel,bf,ttest);
      DVEC dpdb = dbtraj.params().parameters() - refhel.params().parameters();
      for(size_t ipar = 0; ipar < NParams(); ipar++){
	bpgraphs[idir][ipar]->SetPoint(id,dpx[ipar], dpdb[ipar]);
      }
      bgapgraph[idir]->SetPoint(id,dval,(dbtraj.position(ttest)-newbfhel.position(ttest)).R());
    }
    char gtitle[80];
    char gname[80];
    snprintf(gname,80,"db%s",anames[idir].c_str());
    snprintf(gtitle,80,"BFieldMap Change %s",anames[idir].c_str());
    dbcan[idir] = new TCanvas(gname,gtitle,1200,800);
    dbcan[idir]->Divide(3,2);
    for(size_t ipar = 0; ipar < NParams(); ipar++){
      dbcan[idir]->cd(ipar+1);
      bpgraphs[idir][ipar]->Draw("AC*");
    }
    dbcan[idir]->Draw();
    dbcan[idir]->Write();

  }
  TCanvas* dbgcan = new TCanvas("dbgcan","DB Gap",800,800);
  dbgcan->Divide(2,2);
  for(int idir=0;idir<3;++idir){
    dbgcan->cd(idir+1);
    bgapgraph[idir]->Draw("AC*");
  }
  dbgcan->Draw();
  dbgcan->Write();

  lhderiv.Write();
  lhderiv.Close();
  cout << "Return status = " << status << endl;  
  return status;
}

