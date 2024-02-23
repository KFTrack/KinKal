//
// test derivatives of simple trajectory class
//
#include "KinKal/Trajectory/ClosestApproach.hh"

#include <iostream>
#include <cstdio>
#include <iostream>
#include <getopt.h>
#include <typeinfo>

#include "TROOT.h"
#include "TH1F.h"
#include "TFile.h"
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
  printf("Usage: Derivatives  --momentum f --costheta f --azimuth f --particle i --charge i --zorigin f --torigin --delta f --ttest f --By f \n");
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
  KTRAJ reftraj(origin,momv,icharge,bnom);
  //cout << "Reference " << reftraj << endl;
  auto refpos4 = reftraj.position4(ttest);
  //  cout << "origin position " << origin << " test position " << refpos4 << endl;
  auto refmom = reftraj.momentum4(ttest);
  int ndel(50);
  // graphs to compare parameter change
  std::vector<TGraph*> pmomgraphs[3];
  std::vector<TGraph*> pposgraphs[3];
  // graphs to compare momentum change
  TGraph* momgraph[3];
  // gaps
  TGraph* gapgraph[3][3];
  // canvases
  TCanvas* dhcan[3];
  TCanvas* dphcan[3];
  TCanvas* dmomcan[3];
  std::string tfname = KTRAJ::trajName() + "Derivs.root";
  TFile lhderiv(tfname.c_str(),"RECREATE");
  // loop over derivative directions
  double del = 2*delta/(ndel-1);
  double dmin = -delta;
  char gtitle[80];
  char gname[80];
  auto dPdX = reftraj.dPardX(ttest);
  auto dPdM = reftraj.dPardM(ttest);
  auto dXdP = reftraj.dXdPar(ttest);
  auto dMdP = reftraj.dMdPar(ttest);
  // scale of parameter change, for parameter derivative test
  DVEC dpmax;
  for(int idir=0;idir<3;++idir){
    MomBasis::Direction tdir =static_cast<MomBasis::Direction>(idir);
    //    cout << "testing direction " << MomBasis::directionName(tdir) << endl;
    // parameter change
    pmomgraphs[idir] = std::vector<TGraph*>(NParams(),0);
    pposgraphs[idir] = std::vector<TGraph*>(NParams(),0);
    for(size_t ipar = 0; ipar < NParams(); ipar++){
      pmomgraphs[idir][ipar] = new TGraph(ndel);
      string title = KTRAJ::paramName(typename KTRAJ::ParamIndex(ipar));
      title += "#DeltaP;exact;1st derivative";
      pmomgraphs[idir][ipar]->SetTitle(title.c_str());
      pposgraphs[idir][ipar] = new TGraph(ndel);
      title = KTRAJ::paramName(typename KTRAJ::ParamIndex(ipar));
      title += "#DeltaX;exact;1st derivative";
      pposgraphs[idir][ipar]->SetTitle(title.c_str());
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
    DVEC ppder,pmder;

    for(int id=0;id<ndel;++id){
      double dval = dmin + del*id;
      //      cout << "Delta = " << dval << endl;
      // compute 1st order change in parameters
      VEC3 dmomdir = reftraj.direction(ttest,tdir);
      //  compute exact altered params
      VEC3 newmom = refmom.Vect() + dval*dmomdir;
      MOM4 momv(newmom.X(),newmom.Y(),newmom.Z(),pmass);
      KTRAJ xhel(refpos4,momv,icharge,bnom);
      pmder = dPdM*SVEC3(dmomdir.X(), dmomdir.Y(), dmomdir.Z());
      DVEC dvec = reftraj.params().parameters() + dval*pmder;
      // search for max
      for(size_t ipar = 0; ipar < NParams(); ipar++)dpmax[ipar] = std::max(fabs(dvec[ipar]),dpmax[ipar]);
      Parameters pdata(dvec,reftraj.params().covariance());
      KTRAJ dhel(pdata,reftraj);
      // test
      //      auto xpos = xhel.position4(ttest);
      auto dpos = dhel.position4(ttest);
      //      cout << " exa pos " << xpos << endl
      //      << " del pos " << dpos << endl;
      MOM4 dmom = dhel.momentum4(ttest);
      //      cout << "Exact change" << xhel << endl;
      //      cout << "Derivative  " << dhel << endl;
      VEC4 gap = dpos - refpos4;
      // project along 3 directions
      for(int jdir=0;jdir < 3;jdir++){
        MomBasis::Direction tjdir =static_cast<MomBasis::Direction>(jdir);
        VEC3 jmomdir = reftraj.direction(ttest,tjdir);
        gapgraph[idir][jdir]->SetPoint(id,dval,gap.Vect().Dot(jmomdir));
      }
      // parameter diff
      for(size_t ipar = 0; ipar < NParams(); ipar++){
        pmomgraphs[idir][ipar]->SetPoint(id,xhel.paramVal(ipar)-reftraj.paramVal(ipar),dhel.paramVal(ipar)-reftraj.paramVal(ipar));
      }
      // compare momenta after change
      //
      VEC3 dxmom = momv.Vect() - refmom.Vect();
      VEC3 ddmom = dmom.Vect() - refmom.Vect();
      momgraph[idir]->SetPoint(id,dxmom.Dot(dmomdir),ddmom.Dot(dmomdir));
      // now same for position
      auto newpos4 = refpos4 + VEC4(dval*dmomdir.X(),dval*dmomdir.Y(),dval*dmomdir.Z(),0.0);
      KTRAJ xphel(newpos4,refmom,icharge,bnom);
      ppder = dPdX*SVEC3(dmomdir.X(), dmomdir.Y(), dmomdir.Z());
      dvec = reftraj.params().parameters() + dval*ppder;
      for(size_t ipar = 0; ipar < NParams(); ipar++)dpmax[ipar] = std::max(fabs(dvec[ipar]),dpmax[ipar]);
      pdata = Parameters(dvec,reftraj.params().covariance());
      KTRAJ dphel(pdata,reftraj);
      for(size_t ipar = 0; ipar < NParams(); ipar++){
        pposgraphs[idir][ipar]->SetPoint(id,xphel.paramVal(ipar)-reftraj.paramVal(ipar),dphel.paramVal(ipar)-reftraj.paramVal(ipar));
      }
    }
    snprintf(gname,80,"dhMom%s",MomBasis::directionName(tdir).c_str());
    snprintf(gtitle,80,"KTraj Change momentum %s",MomBasis::directionName(tdir).c_str());
    dhcan[idir] = new TCanvas(gname,gtitle,1200,800);
    dhcan[idir]->Divide(3,2);
    TF1* pline = new TF1("pline","[0]+[1]*x");
    for(size_t ipar = 0; ipar < NParams(); ipar++){
      dhcan[idir]->cd(ipar+1);
      // if this is non-trivial, fit
      if(fabs(pmder[ipar])>1e-9){
        pline->SetParameters(0.0,1.0);
        TFitResultPtr pfitr = pmomgraphs[idir][ipar]->Fit(pline,"SQ","AC*");
        pmomgraphs[idir][ipar]->Draw("AC*");
        if(fabs(pfitr->Parameter(0))> 10*delta || fabs(pfitr->Parameter(1)-1.0) > delta){
          cout << "Momentum derivative for parameter "
            << KTRAJ::paramName(typename KTRAJ::ParamIndex(ipar))
            << " in direction " << MomBasis::directionName(tdir)
            << " Out of tolerance : Offset " << pfitr->Parameter(0) << " Slope " << pfitr->Parameter(1) << endl;
          status = 1;
        }
      }
      pmomgraphs[idir][ipar]->Draw("AC*");
    }
    dhcan[idir]->Draw();
    dhcan[idir]->Write();

    snprintf(gname,80,"dhPos%s",MomBasis::directionName(tdir).c_str());
    snprintf(gtitle,80,"KTraj Change position %s",MomBasis::directionName(tdir).c_str());
    dphcan[idir] = new TCanvas(gname,gtitle,1200,800);
    dphcan[idir]->Divide(3,2);
    pline = new TF1("pline","[0]+[1]*x");
    for(size_t ipar = 0; ipar < NParams(); ipar++){
      dphcan[idir]->cd(ipar+1);
      // if this is non-trivial, fit
      if(fabs(ppder[ipar])>1e-9){
        pline->SetParameters(0.0,1.0);
        TFitResultPtr pfitr = pposgraphs[idir][ipar]->Fit(pline,"SQ","AC*");
        pposgraphs[idir][ipar]->Draw("AC*");
        if(fabs(pfitr->Parameter(0))> 10*delta || fabs(pfitr->Parameter(1)-1.0) > delta){
          cout << "Position deriviative for parameter "
            << KTRAJ::paramName(typename KTRAJ::ParamIndex(ipar))
            << " in direction " << MomBasis::directionName(tdir)
            << " Out of tolerance : Offset " << pfitr->Parameter(0) << " Slope " << pfitr->Parameter(1) << endl;
          status = 1;
        }
      }
      pposgraphs[idir][ipar]->Draw("AC*");
    }
    dphcan[idir]->Draw();
    dphcan[idir]->Write();
    //
    snprintf(gname,80,"dm%s",MomBasis::directionName(tdir).c_str());
    snprintf(gtitle,80,"Mom Change %s",MomBasis::directionName(tdir).c_str());
    dmomcan[idir] = new TCanvas(gname,gtitle,800,800);
    dmomcan[idir]->Divide(2,2);
    dmomcan[idir]->cd(1);
    pline->SetParameters(0.0,1.0);
    TFitResultPtr pfitr = momgraph[idir]->Fit(pline,"SQ","AC*");
    momgraph[idir]->Draw("AC*");
    if(fabs(pfitr->Parameter(0))> 10*delta || fabs(pfitr->Parameter(1)-1.0) > 0.01*delta){
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

  // now for parameters
  std::vector<TGraph*> mompgraphs[NParams()];
  std::vector<TGraph*> pospgraphs[NParams()];

  auto refmom3 = reftraj.momentum3(ttest);
  auto refpos3 = reftraj.position3(ttest);
  TCanvas* dparcan[NParams()];
  std::vector<std::string> dirnames = {"X", "Y", "Z"};
  for(size_t ipar = 0; ipar < NParams(); ipar++){
    mompgraphs[ipar] = std::vector<TGraph*>(3,0);
    pospgraphs[ipar] = std::vector<TGraph*>(3,0);
    snprintf(gname,80,"dp%s",KTRAJ::paramName(typename KTRAJ::ParamIndex(ipar)).c_str());
    snprintf(gtitle,80,"Position and Momentum Change WRT %s",KTRAJ::paramName(typename KTRAJ::ParamIndex(ipar)).c_str());
    dparcan[ipar] = new TCanvas(gname,gtitle,1200,800);
    dparcan[ipar]->Divide(3,2);
    for(int idir=0;idir<3;++idir){
      mompgraphs[ipar][idir] = new TGraph(ndel);
      MomBasis::Direction tdir =static_cast<MomBasis::Direction>(idir);
      string title = "Momentum #Delta " + MomBasis::directionName(tdir) + ";exact;1st derivative";
      mompgraphs[ipar][idir]->SetTitle(title.c_str());
      pospgraphs[ipar][idir] = new TGraph(ndel);
      title = "Position #Delta " + MomBasis::directionName(tdir) + ";exact;1st derivative";
      pospgraphs[ipar][idir]->SetTitle(title.c_str());
    }
    for(int id=0;id<ndel;++id){
      double dval = dmin + del*id;
      DVEC newpars = reftraj.params().parameters();
      newpars[ipar] += dpmax[ipar]*dval;
      DVEC dpars = newpars - reftraj.params().parameters();
      Parameters pdata(newpars,reftraj.params().covariance());
      KTRAJ dphel(pdata,reftraj);
      auto dxmom = dphel.momentum3(ttest) - refmom3;
      auto dxpos = dphel.position3(ttest) - refpos3;
      // now derivatives
      SVEC3 dpos = dXdP*dpars;
      SVEC3 dmom = dMdP*dpars;
      VEC3 ddpos(dpos[0], dpos[1], dpos[2]);
      VEC3 ddmom(dmom[0], dmom[1], dmom[2]);
      // project differences along the mom bases
      for(int idir=0;idir<3;++idir){
        MomBasis::Direction tdir =static_cast<MomBasis::Direction>(idir);
        VEC3 jdir = reftraj.direction(ttest,tdir);
        double dxposd = jdir.Dot(dxpos);
        double ddposd = jdir.Dot(ddpos);
        pospgraphs[ipar][idir]->SetPoint(id,dxposd,ddposd);
        double dxmomd = jdir.Dot(dxmom);
        double ddmomd = jdir.Dot(ddmom);
        mompgraphs[ipar][idir]->SetPoint(id,dxmomd,ddmomd);
      }
    }
    TF1* pline = new TF1("pline","[0]+[1]*x");
    for(size_t idir = 0; idir < 3; idir++){
      MomBasis::Direction tdir =static_cast<MomBasis::Direction>(idir);
      VEC3 jdir = reftraj.direction(ttest,tdir);
      SVEC3 jvec(jdir.X(),jdir.Y(),jdir.Z());
      DVEC dp = jvec*dXdP;
      DVEC dm = jvec*dMdP;
      dparcan[ipar]->cd(idir+1);
      // exclude quadratic terms
      double pdiff = pospgraphs[ipar][idir]->GetPointX(ndel-1)-pospgraphs[ipar][idir]->GetPointX(0);
      double pmid = pospgraphs[ipar][idir]->GetPointX(ndel/2-1)-pospgraphs[ipar][idir]->GetPointX(0);
      if(fabs(dp[ipar])>1e-4 && fabs(pdiff)>fabs(pmid)){
        pline->SetParameters(0.0,1.0);
        TFitResultPtr pfitr = pospgraphs[ipar][idir]->Fit(pline,"SQ","AC*");
        if( fabs(pfitr->Parameter(0))> 10*delta || fabs(pfitr->Parameter(1)-1.0) > delta){
          cout << "dXdP for parameter "
            << KTRAJ::paramName(typename KTRAJ::ParamIndex(ipar))
            << " in direction " << MomBasis::directionName(tdir)
            << " Out of tolerance : Offset " << pfitr->Parameter(0) << " Slope " << pfitr->Parameter(1)
            << " derivative  " << dp[ipar] << endl;
          status = 1;
        }
      }
      pospgraphs[ipar][idir]->Draw("AC*");

      dparcan[ipar]->cd(idir+4);
      if(fabs(dm[ipar])>1e-6){
        pline->SetParameters(0.0,1.0);
        TFitResultPtr pfitr = mompgraphs[ipar][idir]->Fit(pline,"SQ","AC*");
        if(fabs(pfitr->Parameter(0))> 10*delta || fabs(pfitr->Parameter(1)-1.0) > delta){
          cout << "dMdP for parameter "
            << KTRAJ::paramName(typename KTRAJ::ParamIndex(ipar))
            << " in direction " << MomBasis::directionName(tdir)
            << " Out of tolerance : Offset " << pfitr->Parameter(0) << " Slope " << pfitr->Parameter(1)
            << " derivative  " << dm[ipar] << endl;
          status = 1;
        }
      }
      mompgraphs[ipar][idir]->Draw("AC*");
    }
    dparcan[ipar]->Draw();
    dparcan[ipar]->Write();
  }

  // test parameter<->phase space translation
  auto dPdS = reftraj.dPardState(ttest);
  auto dSdP = reftraj.dStatedPar(ttest);
  auto ptest = dPdS*dSdP;
  for(size_t irow=0;irow<NParams();irow++) {
    for(size_t icol=0;icol<NParams();icol++) {
      double val(0.0);
      if(irow==icol)val = 1.0;
      if(fabs(ptest(irow,icol) - val) > 1e-9){
        cout <<"Error in parameter derivative test, row col = " << KTRAJ::paramName(typename KTRAJ::ParamIndex(irow))
          << " " << KTRAJ::paramName(typename KTRAJ::ParamIndex(icol))
          << " diff = " << ptest(irow,icol) - val << endl;
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
        cout <<"Error in position derivative test, row col = " << KTRAJ::paramName(typename KTRAJ::ParamIndex(irow))
          << " " << KTRAJ::paramName(typename KTRAJ::ParamIndex(icol))
          << " diff = " << xtest(irow,icol) - val << endl;
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
        cout <<"Error in momentum derivative test"
          << " row " << KTRAJ::paramName(typename KTRAJ::ParamIndex(irow))
          << " col " << KTRAJ::paramName(typename KTRAJ::ParamIndex(icol))
          << " diff = " << mtest(irow,icol) - val << endl;
        status = 1;
      }
    }
  }
  if(status ==1) {
    cout << " mtest" << endl
      << mtest << endl;
    cout << " ptest " << endl
      << ptest << endl;
    cout << " dMdP" << endl
      << dMdP << endl;
    cout << " dPdM" << endl
      << dPdM << endl;
  }

  // test changes due to BFieldMap
  TCanvas* dbcan[3]; // 3 directions
  std::vector<TGraph*> bgraphs[3];
  std::array<VEC3,3> basis;
  basis[0] = reftraj.bnom().Unit();
  basis[1] = reftraj.direction(MomBasis::phidir_);
  basis[2] = VEC3(basis[1].Y(),-basis[1].X(), 0.0); //perp
  std::array<std::string,3> bnames{"BDirection", "PhiDirection", "MomPerpendicular"};
  // gaps
  TGraph* bgapgraph[3];
  auto state = reftraj.stateEstimate(ttest);
  for(size_t idir =0; idir<3;idir++){
    bgraphs[idir] = std::vector<TGraph*>(NParams(),0);
    for(size_t ipar = 0; ipar < NParams(); ipar++){
      bgraphs[idir][ipar] = new TGraph(ndel);
      string title = KTRAJ::paramName(typename KTRAJ::ParamIndex(ipar));
      title += ";exact;1st derivative";
      bgraphs[idir][ipar]->SetTitle(title.c_str());
    }
    bgapgraph[idir] = new TGraph(ndel);
    string title = "Gap for #Delta B in " + bnames[idir] + ";Fractional change;Gap value (mm)";
    bgapgraph[idir]->SetTitle(title.c_str());
    DVEC dpx, dpdb;
    TF1* pline = new TF1("pline","[0]+[1]*x");
    for(int id=0;id<ndel;++id){
      // construct exact helix for this field and the corresponding exact parameter change
      double dval = dmin + del*id;
      VEC3 bf = bnom + basis[idir]*dval/10.0;
      // exact traj given the full state
      KTRAJ newbfhel(state,bf);
      auto newstate = newbfhel.stateEstimate(ttest);
      for(size_t ipar=0;ipar < ParticleState::dimension(); ipar++){
        if(fabs(state.state()[ipar] - newstate.state()[ipar])>1.0e-6) cout << "State vector " << ipar << " doesn't match: original "
          << state.state()[ipar] << " rotated " << newstate.state()[ipar]  << endl;
      }
      dpx = newbfhel.params().parameters() - reftraj.params().parameters();
      // 1st order change trajectory
      KTRAJ dbtraj(reftraj,bf,ttest);
      dpdb = dbtraj.params().parameters() - reftraj.params().parameters();
      for(size_t ipar = 0; ipar < NParams(); ipar++){
        bgraphs[idir][ipar]->SetPoint(id,dpx[ipar], dpdb[ipar]);
      }
      bgapgraph[idir]->SetPoint(id,dval,(dbtraj.position3(ttest)-newbfhel.position3(ttest)).R());
      // BField derivatives

    }
    char gtitle[80];
    char gname[80];
    snprintf(gname,80,"db%s",bnames[idir].c_str());
    snprintf(gtitle,80,"BFieldMap Change %s",bnames[idir].c_str());
    dbcan[idir] = new TCanvas(gname,gtitle,1200,800);
    dbcan[idir]->Divide(3,2);
    for(size_t ipar = 0; ipar < NParams(); ipar++){
      dbcan[idir]->cd(ipar+1);
      pline->SetParameters(0.0,1.0);
      // exclude quadratic terms
      double pdiff = bgraphs[idir][ipar]->GetPointX(ndel-1)-bgraphs[idir][ipar]->GetPointX(0);
      double pmid = bgraphs[idir][ipar]->GetPointX(ndel/2-1)-bgraphs[idir][ipar]->GetPointX(0);
      if(fabs(dpdb[ipar])>1e-9 && fabs(pdiff)>fabs(pmid)){
        TFitResultPtr pfitr = bgraphs[idir][ipar]->Fit(pline,"SQ","AC*");
        if(fabs(pfitr->Parameter(0))> 10*delta || fabs(pfitr->Parameter(1)-1.0) > delta){
          cout << "BField change derivative for parameter "
            << KTRAJ::paramName(typename KTRAJ::ParamIndex(ipar))
            << " in direction " << bnames[idir]
            << " Out of tolerance : Offset " << pfitr->Parameter(0) << " Slope " << pfitr->Parameter(1) << endl;
          status = 1;
        }
      }
      bgraphs[idir][ipar]->Draw("AC*");
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

