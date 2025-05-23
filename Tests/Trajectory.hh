//
// test basic functions of kinematic trajectory class
//
#include "KinKal/Trajectory/SensorLine.hh"
#include "KinKal/Trajectory/ClosestApproach.hh"
#include "KinKal/General/ParticleStateEstimate.hh"
#include "KinKal/General/PhysicalConstants.h"

#include <iostream>
#include <cstdio>
#include <iostream>
#include <getopt.h>

#include "TH1F.h"
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
#include "TProfile.h"
#include "TProfile2D.h"

using namespace KinKal;
using namespace std;

void print_usage() {
  printf("Usage: Trajectory --momentum f --costheta f --azimuth f --particle i --charge i --xorigin f -- yorigin f --zorigin f --torigin --tmin f--tmax f --ltime f --By f --invert i\n");
}

struct MomVec {
  TPolyLine3D* arrow;
  TPolyMarker3D *start, *end;
  MomVec() : arrow(new TPolyLine3D(2)), start(new TPolyMarker3D(1,21)), end(new TPolyMarker3D(1,22)) {}
};

void drawMom(VEC3 const& start, VEC3 const& momvec,int momcolor,MomVec& mom) {
  mom.arrow->SetPoint(0,start.X(),start.Y(),start.Z());
  auto end = start + momvec;
  mom.arrow->SetPoint(1,end.X(),end.Y(),end.Z());
  mom.arrow->SetLineColor(momcolor);
  mom.arrow->Draw();
  mom.start->SetPoint(1,start.X(),start.Y(),start.Z());
  mom.start->SetMarkerColor(momcolor);
  mom.start->Draw();
  mom.end->SetPoint(1,end.X(),end.Y(),end.Z());
  mom.end->SetMarkerColor(momcolor);
  mom.end->Draw();
}

template <class KTRAJ>
int TrajectoryTest(int argc, char **argv,KinKal::DVEC sigmas) {
  int opt;
  double mom(105.0), cost(0.7), phi(0.5);
  double masses[5]={0.511,105.66,139.57, 493.68, 938.0};
  int imass(0), icharge(-1);
  double pmass, ox(0.0), oy(10.0), oz(100.0), ot(0.0);
  double tmin(-5.0), tmax(5.0);
  double ltime(3.0), vprop(0.8), gap(2.0);
  double wlen(1000.0); // length of the wire
  double By(0.0);
  int invert(0);

  static struct option long_options[] = {
    {"momentum",     required_argument, 0, 'm' },
    {"costheta",     required_argument, 0, 'c'  },
    {"azimuth",     required_argument, 0, 'a'  },
    {"particle",     required_argument, 0, 'p'  },
    {"charge",     required_argument, 0, 'q'  },
    {"xorigin",     required_argument, 0, 'x'  },
    {"yorigin",     required_argument, 0, 'y'  },
    {"zorigin",     required_argument, 0, 'z'  },
    {"torigin",     required_argument, 0, 't'  },
    {"tmin",     required_argument, 0, 's'  },
    {"tmax",     required_argument, 0, 'e'  },
    {"ltime",     required_argument, 0, 'l'  },
    {"By",     required_argument, 0, 'b'  },
    {"invert",     required_argument, 0, 'I'  },
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
      case 'x' : ox = atof(optarg);
                 break;
      case 'y' : oy = atof(optarg);
                 break;
      case 'z' : oz = atof(optarg);
                 break;
      case 't' : ot = atof(optarg);
                 break;
      case 's' : tmin = atof(optarg);
                 break;
      case 'e' : tmax = atof(optarg);
                 break;
      case 'l' : ltime = atof(optarg);
                 break;
      case 'b' : By = atof(optarg);
                 break;
      case 'I' : invert = atoi(optarg);
                 break;
      default: print_usage();
               exit(EXIT_FAILURE);
    }
  }

  pmass = masses[imass];

  printf("Testing Kinematic Trajectory %s with momentum = %f, costheta = %f, phi = %f, mass = %f, charge = %i, z = %f, t = %f \n",KTRAJ::trajName().c_str(),mom,cost,phi,pmass,icharge,oz,ot);
  // define the BF (tesla)
  VEC3 bnom(0.0,By,1.0);
  VEC4 origin(ox,oy,oz,ot);
  double sint = sqrt(1.0-cost*cost);
  MOM4 momv(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost,pmass);
  KTRAJ ktraj(origin,momv,icharge,bnom,TimeRange(-10,10));

// test state
  auto state = ktraj.state(3.0);
  KTRAJ straj(state,bnom,ktraj.range());
  for(size_t ipar=0; ipar < NParams(); ipar++){
    if(fabs(ktraj.paramVal(ipar)-straj.paramVal(ipar)) > 1.0e-10){
      std::cout << "Parameter mismatch, par " << ipar << " diff " <<  ktraj.paramVal(ipar) << " " << straj.paramVal(ipar) << std::endl;
      return -1;
    }
  }

  if(invert)ktraj.invertCT();
  auto testmom = ktraj.momentum4(ot);
  //  cout << "KTRAJ with momentum " << momv.Vect() << " position " << origin << " has parameters: " << ktraj << endl;
  //  cout << "origin time position = " << ktraj.position3(ot) << " momentum " << ktraj.momentum3(ot) << " mag " <<  ktraj.momentum(ot) << endl;
  VEC3 tvel, tdir;
  double ttime;
  double tstp = ktraj.range().range()/9;
  for(int istep=0;istep<10;istep++){
    ttime = ktraj.range().begin() + istep*tstp;
    tvel = ktraj.velocity(ttime);
    tdir = ktraj.direction(ttime);
    testmom = ktraj.momentum4(ttime);
    //    cout << "velocity " << tvel << " direction " << tdir << " momentum " << testmom.R() << endl;
    //    cout << "momentum beta =" << testmom.Beta() << " KTRAJ beta = " << ktraj.beta() << " momentum gamma  = " << testmom.Gamma() <<
    //      " KTRAJ gamma = " << ktraj.gamma() << " scalar mom " << ktraj.momentum(ot) << endl;
  }
  VEC3 mdir = ktraj.direction(ot);
  // create the helix at tmin and tmax
  auto tmom = ktraj.momentum4(tmax);
  auto tpos = ktraj.position4(tmax);
  KTRAJ ktrajmax(tpos,tmom,icharge,bnom);
  tmom = ktraj.momentum4(tmin);
  tpos = ktraj.position4(tmin);
  KTRAJ ktrajmin(tpos,tmom,icharge,bnom);

  //  cout << "KTRAJ at tmax has parameters : " << ktrajmax << endl;
  //  cout << "KTRAJ at tmin has parameters : " << ktrajmin << endl;

  // now graph this as a polyline over the specified time range.
  double tstep = 0.1; // nanoseconds
  double trange = tmax-tmin;
  int nsteps = (int)rint(trange/tstep);
  // create Canvase
  TCanvas* hcan = new TCanvas("hcan","Helix",1000,1000);
  //TPolyLine to graph the result
  TPolyLine3D* hel = new TPolyLine3D(nsteps+1);
  for(int istep=0;istep<nsteps+1;++istep){
    // compute the position from the time
    VEC4 hpos = ktraj.position4(tmin + tstep*istep);
    // add these positions to the TPolyLine3D
    hel->SetPoint(istep, hpos.X(), hpos.Y(), hpos.Z());
  }
  // draw the helix
  if(icharge > 0)
    hel->SetLineColor(kBlue);
  else
    hel->SetLineColor(kRed);
  hel->Draw();
  // inversion test
  TPolyLine3D* ihel = new TPolyLine3D(nsteps+1);
  auto iktraj = ktraj;
  iktraj.invertCT();
  for(int istep=0;istep<nsteps+1;++istep){
    // compute the position from the time
    VEC4 hpos = ktraj.position4(tmin + tstep*istep);
    // add these positions to the TPolyLine3D
    ihel->SetPoint(istep, hpos.X(), hpos.Y(), hpos.Z());
  }
  // draw the helix
  ihel->SetLineColor(kBlue);
  ihel->SetLineStyle(kDashDotted);
  ihel->Draw();
  // now draw momentum vectors at reference, start and end
  MomVec imstart,imend,imref;
  auto imompos = iktraj.position3(ot);
  mdir =iktraj.direction(ot);
  VEC3 imomvec =mom*mdir;
  drawMom(imompos,imomvec,kRed,imref);
  //
  imompos  = ktraj.position3(tmin);
  mdir = ktraj.direction(tmin);
  imomvec =mom*mdir;
  drawMom(imompos,imomvec,kBlue,imstart);
  //
  imompos = ktraj.position3(tmax);
  mdir = ktraj.direction(tmax);
  imomvec =mom*mdir;
  drawMom(imompos,imomvec,kGreen,imend);
  //

  // draw the origin and axes
  TAxis3D* rulers = new TAxis3D();
  rulers->GetXaxis()->SetAxisColor(kBlue);
  rulers->GetXaxis()->SetLabelColor(kBlue);
  rulers->GetYaxis()->SetAxisColor(kCyan);
  rulers->GetYaxis()->SetLabelColor(kCyan);
  rulers->GetZaxis()->SetAxisColor(kOrange);
  rulers->GetZaxis()->SetLabelColor(kOrange);
  rulers->Draw();

  // now draw momentum vectors at reference, start and end
  MomVec mstart,mend,mref;
  VEC3 mompos = ktraj.position3(ot);
  mdir = ktraj.direction(ot);
  VEC3 momvec =mom*mdir;
  drawMom(mompos,momvec,kBlack,mref);
  //
  mompos = ktraj.position3(tmin);
  mdir = ktraj.direction(tmin);
  momvec =mom*mdir;
  drawMom(mompos,momvec,kBlue,mstart);
  //
  mompos = ktraj.position3(tmax);
  mdir = ktraj.direction(tmax);
  momvec =mom*mdir;
  drawMom(mompos,momvec,kGreen,mend);
  //
  TLegend* leg = new TLegend(0.8,0.8,1.0,1.0);
  char title[80];
  snprintf(title,80,"KTRAJ, q=%1i, mom=%3.1g MeV/c",icharge,mom);
  leg->AddEntry(hel,title,"L");
  snprintf(title,80,"Ref. Momentum, t=%4.2g ns",ot);
  leg->AddEntry(mref.arrow,title,"L");
  snprintf(title,80,"Inverted Momentum, t=%4.2g ns",ot);
  leg->AddEntry(imref.arrow,title,"L");
  snprintf(title,80,"Start Momentum, t=%4.2g ns",ot+tmin);
  leg->AddEntry(mstart.arrow,title,"L");
  snprintf(title,80,"End Momentum, t=%4.2g ns",ot+tmax);
  leg->AddEntry(mend.arrow,title,"L");
  leg->Draw();

  // create a Line near this helix, and draw it and the ClosestApproach vector
  VEC3 pos = ktraj.position3(ltime);
  VEC3 dir = ktraj.direction(ltime);
  // rotate the direction
  double lhphi = atan2(dir.Y(),dir.X());
  double pphi = lhphi + M_PI/2.0;
  VEC3 pdir(cos(pphi),sin(pphi),0.0);
  double pspeed = CLHEP::c_light*vprop; // vprop is relative to c
  VEC3 pvel = pdir*pspeed;
  // shift the position
  VEC3 perpdir(-sin(phi),cos(phi),0.0);
  VEC3 ppos = pos + gap*perpdir;
  SensorLine tline(ppos, ltime, pvel, wlen);
  // find ClosestApproach
  CAHint hint(ltime,ltime);
  ClosestApproach<KTRAJ,SensorLine> tp(ktraj,tline,hint, 1e-6);
  //  cout << "ClosestApproach status " << tp.statusName() << " doca " << tp.doca() << " dt " << tp.deltaT() << endl;
  if(tp.status() == ClosestApproachData::converged) {
    // draw the line and ClosestApproach
    TPolyLine3D* line = new TPolyLine3D(2);
    auto plow = tline.start();
    auto phigh = tline.end();
    line->SetPoint(0,plow.X(),plow.Y(), plow.Z());
    line->SetPoint(1,phigh.X(),phigh.Y(), phigh.Z());
    line->SetLineColor(kOrange);
    line->Draw();
    TPolyLine3D* poca = new TPolyLine3D(2);
    poca->SetPoint(0,tp.particlePoca().X() ,tp.particlePoca().Y() ,tp.particlePoca().Z());
    poca->SetPoint(1,tp.sensorPoca().X() ,tp.sensorPoca().Y() ,tp.sensorPoca().Z());
    poca->SetLineColor(kBlack);
    poca->Draw();
  } else {
    cout << "ClosestApproach failed " << endl;
    return -1;
  }
  // test particle state back-and-forth
  // set the covariance first
  for(size_t ipar=0; ipar < NParams(); ipar++){
    ktraj.params().covariance()[ipar][ipar] = sigmas[ipar]*sigmas[ipar];
  }
  ParticleStateEstimate pmeas = ktraj.stateEstimate(ltime);
  KTRAJ newktraj(pmeas,ktraj.bnom());
  for(size_t ipar=0;ipar < NParams();ipar++){
    if(fabs(ktraj.paramVal(ipar)-newktraj.paramVal(ipar)) > 1e-9){
      cout << "Parameter check failed par " << ipar << endl;
      return -2;
    }
    if(fabs(ktraj.paramVar(ipar)-newktraj.paramVar(ipar)) > 1e-9){
      cout << "Parameter covariance check failed par " << ipar << endl;
      return -2;
    }
  }
  // test position and momentum variance
  double pmvar = pmeas.momentumVariance();
  double pphivar = pmeas.positionVariance(MomBasis::phidir_);
  double pperpvar = pmeas.positionVariance(MomBasis::perpdir_);

  double tmvar = ktraj.momentumVariance(ltime);
  double tphivar = ktraj.positionVariance(ltime,MomBasis::phidir_);
  double tperpvar = ktraj.positionVariance(ltime,MomBasis::perpdir_);

  auto momdir = ktraj.direction(ltime);
  auto udir = ktraj.direction(ltime,MomBasis::perpdir_);
  Plane plane(momdir,udir,ppos);
  auto pcov = ktraj.planeCovariance(ltime,plane);
//  cout << "Plane covariance " << pcov << endl;

  if(fabs(pcov(0,0) - pperpvar) > 1e-9 ||
      fabs(pcov(1,1) - pphivar) > 1e-9) {
    cout << "Position covariance check failed" << endl;
    return -2;
  }

  // test acceleration
  auto acc = ktraj.acceleration(ltime);
  auto vel = ktraj.velocity(ltime);
  if(acc.Dot(vel) > 1e-9 || acc.Dot(ktraj.bnom()) > 1e-9){
    cout << "Acceleration check failed " << endl;
    return -2;
  }

  if(fabs(pmvar-tmvar)>1e-9 ) {
    cout << "Momentum covariance check failed " << endl;
    return -2;
  }

  if(fabs(pphivar-tphivar) > 1e-9 || fabs(pperpvar-tperpvar) > 1e-9){
    cout << "Position covariance check failed " << endl;
    return -2;
  }

  std::string tfname = KTRAJ::trajName() + ".root";
  cout << "Saving canvas to " << title << endl;
  hcan->SaveAs(tfname.c_str());
  cout <<"Exiting with status " << EXIT_SUCCESS << endl;
  exit(EXIT_SUCCESS);
}

