//
//  Class definitions and tests for a kinematic 6-parameter track fit
//  Original author: David Brown (LBNL), Dec. 2018
//
#include "TH1F.h"
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
#include "TMath.h"
#include "Math/SVector.h"
#include "Math/SMatrix.h"
//#include "Math/PxPyPzE4D.h"
#include <math.h>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace ROOT::Math;

// global values: these should be obtained from services 
enum effect {momfrac=0,theta1,theta2};
enum pars {rad_=0, lam_=1, cx_=2, cy_=3, phi0_=4, t0_=5, _npars=6};
vector<string> parnames = {"Radius", "Lambda", "Cx", "Cy",  "Phi0", "T0"};
vector<string> parunits = {"(mm)", "(mm)", "(mm)", "(mm)",  "", "(nsec)"};
//  rad_; // transverse radius
//  lam_; // longiduinal wavelength
//  cx_, cy_; // cylinder center transverse coordinates
//  phi0_; // phi angle at z=0 plane
//  t0_; // time helix passes z=0 plane
double c=299.792;
double twopi = 6.283185;
double pi = twopi/2.0;
// set B to 1.0
double B=1.0;
// drift velocity 
double vdrift(0.06); // drift velocity
double timez = 1600.0; // location of time hit

double minDOCA = 0.0;  // don't use drift for DOCA below this
double rmin(350.0), rmax(700.0); // radial limits
TRandom* gRandom = new TRandom3(324913);

typedef SVector<double,6> PVec;
typedef SMatrix<double,6,6,MatRepSym<double,6> > PMat;
typedef SMatrix<double,6,1> PDer;

//struct FourV {
//  double _x, _y, _z, _t;
//  FourV(): _x(0), _y(0), _z(0), _t(0) {}
//  FourV(double x, double y, double z, double t): _x(x), _y(y), _z(z), _t(t) {}
//  double dot3(FourV& other) {
//    return _x*other._x + _y*other._y + _z*other._z; }
//  double dot4(FourV& other) {
//    return dot3(other) - _t*other._t; }
//  double mag3() { return sqrt(dot3(*this)); }
//  double mag4() { return sqrt(dot4(*this)); }
//  FourV& normalize() { double norm = 1.0/mag3(); _x *= norm; _y *= norm; _z *= norm; return *this; }
//  FourV& operator += (FourV const& other) { _x += other._x; _y += other._y; _z += other._z; _t += other._t; return *this; }
//  FourV& operator -= (FourV const& other) { _x -= other._x; _y -= other._y; _z -= other._z; _t -= other._t; return *this; }
//  friend FourV operator + (FourV left, FourV const& right) { return left += right; }
//  friend FourV operator - (FourV left, FourV const& right) { return left -= right; }
//  // set energy given mass in MeV/c^2
//  void SetEnergy(double mass) { _t = sqrt(dot3(*this) + mass*mass); }
//};
  typedef ROOT::Math::LorentzVector<ThreeV > FourV;

ostream& operator << (ostream& os, FourV const& vec) {
  os << vec._x << " " 
  << vec._y << " " 
  << vec._z << " " 
  << vec._t << endl;
  return os;
}

struct HelixPars {
  PVec _pars; // parameter vector
  // the following is not a true parameter, but it defines the connection between time as a parametric variable
  // and the geometric helix description, so must be included
  double _mbar;  // reduced mass in units of mm
  double _mass;  // in units of MeV
  double pbar() const { return  sqrt(_pars[rad_]*_pars[rad_] + _pars[lam_]*_pars[lam_] ); } // momentum in mm
  double ebar() const { return  sqrt(_pars[rad_]*_pars[rad_] + _pars[lam_]*_pars[lam_] + _mbar*_mbar); } // energy in mm
// angular rotation frequency
  double omega() const { return copysign(c,_mbar)/ebar(); } // rotational velocity, sign set by magnetic force 
  double beta() const { return pbar()/ebar(); }
  double phi(double t) const { return omega()*(t - _pars[t0_]) + _pars[phi0_]; }
  double time(double zpos) const { return _pars[t0_] + zpos/(omega()*_pars[lam_]); }
  void invert() {
    _mbar *= -1.0;
    _pars[t0_] *= -1.0;
  }
};
// drift hit along a helix; generated from a particular helix
struct DHit {
  FourV _hpos; // helix position and time
  FourV _wpos; // wire position and time
  double _weta; // hit wire axis azimumth angle
  double _rdrift; // hit true drift radius, signed by ambiguity
  double _t;  // hit measurement time; depends only on generation parameters
  double _sigt; // hit time smearing sigma
};
// timing measurement hit
struct THit {
  FourV _hpos; // helix position and time
  double _t; // hit measurement time
  double _sigt; // hit time smearing sigma
  PVec _dtoca; // time derivatives
};
  
// DOCA information: this of for a particular helix and hit pair
struct HDOCA {
  int _ambig; // reco ambig
  double _doca; // DOCA
  double _toca; // time of closest approach = propagation time + drift time
  PVec _ddoca; // DOCA derivatives
  PVec _dtoca; // time derivatives
};

void CrossProduct(FourV const& a, FourV const& b, FourV& c) {
  c._t = a._t;
  c._x = a._y*b._z - a._z*b._y;
  c._y = a._z*b._x - a._x*b._z;
  c._z = a._x*b._y - a._y*b._x;
}

// set the space components of the vector given the time
void HelixPos(HelixPars const& pars, FourV& pos) {
// compute rotational frequency
  double omega = pars.omega();
// relative time
  double dt = pos._t - pars_[t0_];
// compute azimuthal angle
  double phi = omega*dt + pars_[phi0_];
// now compute position
  pos._x = pars_[cx_] + pars_[rad_]*sin(phi);
  pos._y = pars_[cy_] - pars_[rad_]*cos(phi);
  pos._z = omega*pars_[lam_]*dt;
}


void HelixVel(HelixPars const& pars, FourV& vel){
// compute rotational frequency
  double omega = pars.omega();
// relative time
  double dt = vel._t - pars_[t0_];
// compute azimuthal angle
  double phi = omega*dt + pars_[phi0_];
// now compute position
  vel._x = omega*pars_[rad_]*cos(phi);
  vel._y = omega*pars_[rad_]*sin(phi);
  vel._z = omega*pars_[lam_];
}

void HelixMom(HelixPars const& pars, double time, FourV& mom){
// compute velocity vector
  double phi = pars.phi(time);
  double factor = pars._mass/pars._mbar;
  mom._x = factor * pars_[rad_] * cos(phi);
  mom._y = factor * pars_[rad_] * sin(phi);
  mom._z = factor * pars_[lam_];
  mom.SetEnergy(pars._mass);
}

// convert 4-momentum positions and momentum, plus particle charge, into a helix.
// Z is assumed to be along Z
void HelixFromMom(FourV pos, FourV mom, double charge, double Bz, HelixPars& pars){
// speed of light in mm/msec
// compute some simple useful parameters
  double pt = sqrt(mom._x*mom._x+mom._y*mom._y);
  double phibar = atan2(mom._y,mom._x);
// translation factor from MeV/c to curvature radius in mm
  double momToRad = 1000.0/(charge*Bz*c);
  // mass in units of mev/c^2
  pars._mass = sqrt(mom._t*mom._t - mom._x*mom._x - mom._y*mom._y - mom._z*mom._z);
// reduced mass
  pars._mbar = -pars._mass*momToRad;
// transverse radius of the helix
  pars_[rad_] = -pt*momToRad;
//tan dip
  pars_[lam_] = -mom._z*momToRad;
  // time at z=0
  pars_[t0_] = pos._t - pos._z/(pars.omega()*pars_[lam_]);
// compute winding that miminimizes z1
  double nwind = rint((pos._z/(pars_[lam_]) - phibar)/twopi);
//  cout << "winding number = " << nwind << endl;
// azimuth at z=0
  pars_[phi0_] = phibar - pars.omega()*(pos._t-pars_[t0_]) + twopi*nwind;
// circle center
  pars_[cx_] = pos._x + mom._y*momToRad;
  pars_[cy_] = pos._y - mom._x*momToRad;
}

void HelixChangeDir(HelixPars const& pars, effect eff,FourV& unit) {
  double phi = pars.omega()*(unit._t - pars_[t0_]) + pars_[phi0_];
  double pbar = pars.pbar();
  switch ( eff ) {
    case theta1:
      unit._x = pars_[lam_]*cos(phi)/pbar;
      unit._y =	pars_[lam_]*sin(phi)/pbar;
      unit._z = pars_[rad_]/pbar;
      break;
    case theta2:
      unit._x = sin(phi);
      unit._y = -cos(phi);
      unit._z = 0.0;
      break;
    case momfrac:
      unit._x = pars_[rad_]*cos(phi)/pbar;
      unit._y = pars_[rad_]*sin(phi)/pbar;
      unit._z = -pars_[lam_]/pbar;
      break;
    default:
      cout << "Error: unknown effect!" << eff << endl;
  }
}

void HelixChange(HelixPars const& pars,HelixPars& newpars, double time, double delta, effect eff) {
// initialize to no effect
  newpars = pars;
// compute some useful quantities
  double beta = pars.beta();
  double omega = pars.omega();
  double pbar = pars.pbar();
  double phi = omega*(time - pars_[t0_]) + pars_[phi0_];
// cases
  switch ( eff ) {
    case theta1:
    // polar bending: only momentum and position are unchanged
      newpars_[rad_] += delta*pars_[lam_];
      newpars_[lam_] += -delta*pars_[rad_];
      newpars_[t0_] += delta*(time-pars_[t0_])/pars_[lam_];
      newpars_[phi0_] += delta*omega*(time-pars_[t0_])/pars_[lam_];
      newpars_[cx_] += -delta*pars_[lam_]*sin(phi);
      newpars_[cy_] += delta*pars_[lam_]*cos(phi);
    break;
    case theta2:
    // Azimuthal bending: R, Lambda, t0 are unchanged
      newpars_[phi0_] += copysign(delta,delta*omega)*pbar/pars_[rad_];
      newpars_[cx_] += -copysign(delta,delta*omega)*pbar*cos(phi);
      newpars_[cy_] += -copysign(delta,delta*omega)*pbar*sin(phi);
    break;
    case momfrac:
    // fractional momentum change: position and direction are unchanged
      newpars_[lam_] += delta*pars_[lam_];
      newpars_[rad_] += delta*pars_[rad_];
      newpars_[t0_] += delta*(time-pars_[t0_])*(1.0-beta*beta);
      newpars_[phi0_] += delta*omega*(time-pars_[t0_]);
      newpars_[cx_] += -delta*pars_[rad_]*sin(phi);
      newpars_[cy_] += +delta*pars_[rad_]*cos(phi);
    break;
    default:
      cout << "Error: unknown effect!" << eff << endl;
  }
}

void TestHelix( double charge, double x, double y, double z, double cost, double momphi, double momval=100.0, double mass=0.5,double tmin=-8, double tmax=8) {
// position and momentum vectors
  FourV pos;
  pos._x = x;
  pos._y = y;
  pos._z = z;
  pos._t = 0.5*tmin+0.5*tmax;
  FourV mom;
  double sint = sqrt(1-cost*cost);
  mom._x = momval*sint*cos(momphi);
  mom._y = momval*sint*sin(momphi);
  mom._z = momval*cost;
  mom.SetEnergy(mass);
// construct the helix from these
  HelixPars pars;
  HelixFromMom(pos,mom,charge,B,pars);
  cout << "Helix parameter:" << endl
  << "Transverse Radius = " << pars_[rad_] << endl
  << "Lambda = " << pars_[lam_] << endl
  << "Time-0 = " << pars_[t0_] << endl
  << "phi-0 = " << pars_[phi0_] << endl
  << "Center x = " << pars_[cx_] << endl
  << "Center y = " << pars_[cy_] << endl
  << "mass = " << pars._mass << endl
  << "Reduced mass = " << pars._mbar << endl;

  FourV testmom;
  HelixMom(pars,pos._t,testmom);
  cout << "Original mom " << mom << endl << "Helix mom " << testmom << endl;
  HelixPars invpars = pars;
  invpars.invert();
  // test velocity
  FourV vel;
  vel._t = pos._t;
  HelixVel(pars,vel);
  double dot = testmom.dot3(vel)/c;
  cout << "velocity dot mom = " << dot << endl;

  double velmag = sqrt(vel._x*vel._x + vel._y*vel._y + vel._z*vel._z);
  double beta = velmag/c;
  double gamma = 1.0/sqrt(1.0 - beta*beta);
  double energy = sqrt(momval*momval + mass*mass);

  cout << "beta = " << beta << " gamma = " << gamma << " E/m = " << energy/mass << " P/E = " << momval/energy << endl;

// now graph this as a polyline over the specified time range.
  double tstep = 0.1; // nanoseconds
  double trange = tmax-tmin;
  int nsteps = (int)rint(trange/tstep);
// create Canvase
  TCanvas* hcan = new TCanvas("hcan","Helix",1000,1000);
//TPolyLine to graph the result
  TPolyLine3D* hel = new TPolyLine3D(nsteps+1);
  TPolyLine3D* invhel = new TPolyLine3D(nsteps+1);
  FourV hpos;
  for(int istep=0;istep<nsteps+1;++istep){
  // compute the position from the time
    hpos._t = tmin + tstep*istep;
    HelixPos(pars,hpos);
    // add these positions to the TPolyLine3D
    hel->SetPoint(istep, hpos._x, hpos._y, hpos._z);
    HelixPos(invpars,hpos);
    invhel->SetPoint(istep, hpos._x, hpos._y, hpos._z);
  }
  // draw the helix
  if(charge > 0)
    hel->SetLineColor(kBlue);
  else
    hel->SetLineColor(kRed);
  hel->Draw();
//  invhel->SetLineColor(kCyan);
//  invhel->Draw();
  // draw the origin and axes
  TAxis3D* rulers = new TAxis3D();
  rulers->GetXaxis()->SetAxisColor(kBlue);
  rulers->GetXaxis()->SetLabelColor(kBlue);
  rulers->GetYaxis()->SetAxisColor(kCyan);
  rulers->GetYaxis()->SetLabelColor(kCyan);
  rulers->GetZaxis()->SetAxisColor(kOrange);
  rulers->GetZaxis()->SetLabelColor(kOrange);
//   TAxis3D::ToggleRulers();
//   TAxis3D::ToggleRulers();
//  cout <<"axis = " << rulers << endl;
  rulers->Draw();

//  TEveArrow* start = new TEveArrow(pos._x,pos._y,pos._z,0.0,0.0,0.0);
//  start->Draw();

//  TVector3* start =new TVector3(pos._x,pos._y,pos._z);
//  start->Draw();
  TPolyLine3D* refmom = new TPolyLine3D(2);
  refmom->SetPoint(0,pos._x,pos._y,pos._z);
  double mommag = sqrt(mom._x*mom._x + mom._y*mom._y + mom._z*mom._z);
  double momscale = fabs(pars_[rad_])/mommag;
  refmom->SetPoint(1,pos._x + mom._x*momscale,pos._y + mom._y*momscale ,pos._z + mom._z*momscale);
  int momcolor;
  if(charge>0.0)
    momcolor = kRed;
  else
    momcolor = kBlack;
  refmom->SetLineColor(momcolor);
  refmom->Draw();

  TPolyMarker3D* refp = new TPolyMarker3D(1,24);
  refp->SetMarkerColor(momcolor);
  refp->SetPoint(0,pos._x,pos._y,pos._z);
  refp->Draw();
  TPolyMarker3D* refmomp = new TPolyMarker3D(1,22);
  refmomp->SetPoint(0,pos._x + mom._x*momscale,pos._y + mom._y*momscale ,pos._z + mom._z*momscale);
  refmomp->SetMarkerColor(momcolor);
  refmomp->Draw();

  TPolyMarker3D* helixp = new TPolyMarker3D(1,3);
  helixp->SetMarkerColor(kGreen);
  hpos._t = pos._t;
  HelixPos(pars,hpos);
  helixp->SetPoint(0,hpos._x,hpos._y,hpos._z);
  helixp->Draw();

  TPolyMarker3D* testp = new TPolyMarker3D(1,25);
  testp->SetMarkerColor(kOrange);
  testp->SetPoint(0,hpos._x,hpos._y,hpos._z+2*M_PI*pars_[lam_]);
  testp->Draw();

  TPolyMarker3D* startp = new TPolyMarker3D(1,21);
  startp->SetMarkerColor(kBlue);
  hpos._t = tmin;
  HelixPos(pars,hpos);
  startp->SetPoint(0,hpos._x,hpos._y,hpos._z);
  startp->Draw();

  TPolyMarker3D* endp = new TPolyMarker3D(1,22);
  endp->SetMarkerColor(kBlue);
  hpos._t = tmax;
  HelixPos(pars,hpos);
  endp->SetPoint(0,hpos._x,hpos._y,hpos._z);
  endp->Draw();

  TLegend* leg = new TLegend(0.8,0.8,1.0,1.0);
  char title[80];
  snprintf(title,80,"Helix, mass=%3.1g MeV/c^{2}, q=%1.1f",mass,charge);
  leg->AddEntry(hel,title,"L");
  snprintf(title,80,"Initial Momentum =%3.1g MeV/c",momval);
  leg->AddEntry(refmom,title,"L");
  leg->AddEntry(refp,"Initial Position","P");
  snprintf(title,80,"Helix, t=%4.2g ns",pos._t);
  leg->AddEntry(helixp,title,"P");
  snprintf(title,80,"Helix, t=%4.2g ns",pos._t+tmin);
  leg->AddEntry(startp,title,"P");
  snprintf(title,80,"Helix, t=%4.2g ns",pos._t+tmax);
  leg->AddEntry(endp,title,"P");
  leg->Draw();

}

void TestHelixDerivs( double dmin, double dmax, double time, effect ieff, double charge=1.0, double momval=100.0, double momphi=1.5, double cost=0.7, double mass=0.5,double t0=0.0){
// position and momentum vectors
  FourV pos;
  pos._x = 0.0;
  pos._y = 0.0;
  pos._z = 0.0;
  pos._t = t0;
  FourV mom;
  double sint = sqrt(1-cost*cost);
  mom._x = momval*sint*cos(momphi);
  mom._y = momval*sint*sin(momphi);
  mom._z = momval*cost;
  mom.SetEnergy(mass);
// construct original helix from these
  HelixPars pars;
  HelixFromMom(pos,mom,charge,B,pars);
// position at time of change
  FourV testpos;
  testpos._t = time;
  HelixPos(pars,testpos);
  FourV testmom;
  HelixMom(pars,time,testmom);
  cout << "initial momentum = " << mom._x << " " << mom._y << " " << mom._z << endl;
  cout << "derived momentum = " << testmom._x << " " << testmom._y << " " << testmom._z << endl;
// check
  HelixPars testpars;
  HelixFromMom(testpos,testmom,charge,B,testpars);
  cout << "Delta radius = " << testpars_[rad_] << " " << pars_[rad_] <<
   " Delta lambda = " << testpars_[lam_] << " " << pars_[lam_] <<
   " Delta t0 = " << testpars_[t0_] << " " << pars_[t0_] <<
   " Delta phi0 = " << testpars_[phi0_] << " " << pars_[phi0_] <<
   " Delta Cx = " << testpars_[cx_] << " " << pars_[cx_] <<
   " Delta Cy = " << testpars_[cy_] << " " << pars_[cy_] << endl;

// momentum direction at time of change
  double testphi = atan2(testmom._y,testmom._x);
  double testcost = testmom._z/sqrt(testmom._x*testmom._x + testmom._y*testmom._y + testmom._z*testmom._z);
  cout << "old phi = " << momphi << " new phi = " << testphi << endl;
  cout << "old cost = " << cost << " new cost = " << testcost << endl;
// make graphs to compare each parameter
  int ndel(100);
  TGraph* radgraph = new TGraph(ndel+1);
  radgraph->SetTitle("Radius;exact;1st derivative");
  TGraph* lambdagraph = new TGraph(ndel+1);
  lambdagraph->SetTitle("Lambda;exact;1st derivative");
  TGraph* t0graph = new TGraph(ndel+1);
  t0graph->SetTitle("t_{0};exact;1st derivative");
  TGraph* phi0graph = new TGraph(ndel+1);
  phi0graph->SetTitle("phi_{0};exact;1st derivative");
  TGraph* cxgraph = new TGraph(ndel+1);
  cxgraph->SetTitle("C_{x};exact;1st derivative");
  TGraph* cygraph = new TGraph(ndel+1);
  cygraph->SetTitle("C_{y};exact;1st derivative");
  // graphs to compare momentum change
  TGraph* mom0graph = new TGraph(ndel+1);
  mom0graph->SetTitle("Momentum Direction;exact;1st derivative");
  TGraph* mom1graph = new TGraph(ndel+1);
  mom1graph->SetTitle("Theta Direction;exact;1st derivative");
  TGraph* mom2graph = new TGraph(ndel+1);
  mom2graph->SetTitle("Phi Direction;exact;1st derivative");

  TGraph* gapgraph = new TGraph(ndel+1);
  gapgraph->SetTitle("Gap;change;gap value (mm)");

// scan range of change
  double del = (dmax-dmin)/ndel;
  for(int id=0;id<=ndel+1;++id){
    double delta = dmin + del*id; 
    //  compute exact altered params
    HelixPars xpars;
    FourV newmom;
    if(ieff == theta1){
      double newcost = cos(acos(cost) + delta);
      double newsint = sqrt(1.0-newcost*newcost);
      newmom._x = momval*newsint*cos(testphi);
      newmom._y = momval*newsint*sin(testphi);
      newmom._z = momval*newcost;
      newmom.SetEnergy(mass);
      HelixFromMom(testpos,newmom,charge,B,xpars);
    } else if(ieff == theta2){
      double dphi = delta/sint;
      double newphi = testphi + dphi;
      newmom._x = momval*sint*cos(newphi);
      newmom._y = momval*sint*sin(newphi);
      newmom._z = momval*cost;
      newmom.SetEnergy(mass);
      HelixFromMom(testpos,newmom,charge,B,xpars);
    } else if(ieff == momfrac) {
      double newmomval = momval*(1.0+delta);
      newmom._x = newmomval*sint*cos(testphi);
      newmom._y = newmomval*sint*sin(testphi);
      newmom._z = newmomval*cost;
      newmom.SetEnergy(mass);
      HelixFromMom(testpos,newmom,charge,B,xpars);
    }
    // now, also compute 1st order change in parameters
    HelixPars d1pars;
    HelixChange(testpars,d1pars, testpos._t, delta, ieff);
    FourV d1pos;
    d1pos._t = testpos._t;
    HelixPos(d1pars,d1pos);
    FourV gap = d1pos - testpos;
    gapgraph->SetPoint(id,delta,gap.mag3());
    FourV d1mom;
    HelixMom(d1pars,testpos._t,d1mom);
    radgraph->SetPoint(id,xpars_[rad_],d1pars_[rad_]);
    lambdagraph->SetPoint(id,xpars_[lam_],d1pars_[lam_]);
    t0graph->SetPoint(id,xpars_[t0_],d1pars_[t0_]);
    phi0graph->SetPoint(id,xpars_[phi0_],d1pars_[phi0_]);
    cxgraph->SetPoint(id,xpars_[cx_],d1pars_[cx_]);
    cygraph->SetPoint(id,xpars_[cy_],d1pars_[cy_]);
    // compare momenta after change
    FourV dxmom = newmom - testmom;
    FourV dd1mom = d1mom - testmom;
    FourV changedir;
    changedir._t = time;
    HelixChangeDir(testpars,momfrac,changedir);
    mom0graph->SetPoint(id,dxmom.dot3(changedir),dd1mom.dot3(changedir));
    HelixChangeDir(pars, theta1, changedir );
    mom1graph->SetPoint(id,dxmom.dot3(changedir),dd1mom.dot3(changedir));
    HelixChangeDir(pars, theta2, changedir );
    mom2graph->SetPoint(id,dxmom.dot3(changedir),dd1mom.dot3(changedir));
  }
  // draw comparisons
  TCanvas* dhcan = new TCanvas("dhcan","Helix Change",1200,800);
  dhcan->Divide(3,2);
  dhcan->cd(1);
  radgraph->Draw("AC*");
  dhcan->cd(2);
  lambdagraph->Draw("AC*");
  dhcan->cd(3);
  t0graph->Draw("AC*");
  dhcan->cd(4);
  phi0graph->Draw("AC*");
  dhcan->cd(5);
  cxgraph->Draw("AC*");
  dhcan->cd(6);
  cygraph->Draw("AC*");

  TCanvas* dmomcan = new TCanvas("dmomcan","Mom Change",800,800);
  dmomcan->Divide(2,2);
  dmomcan->cd(1);
  mom0graph->Draw("AC*");
  dmomcan->cd(2);
  mom1graph->Draw("AC*");
  dmomcan->cd(3);
  mom2graph->Draw("AC*");
  dmomcan->cd(4);
  gapgraph->Draw("AC*");
}

void TestHelixAccuracy(double mommin, double mommax, unsigned nbins, double mass, double xfrac,double costmin, double costmax, double xrange, unsigned nsamples, double rmax, double zmax) {
// set B to 1.0
  double charge = 1.0;
// random generator
// diagnostics
  double maxgap(100.0);
  TProfile2D* mfgapP = new TProfile2D("mfgapP","Energy straggling gap;Momentum (MeV/c);cos(#theta);gap (mm)",nbins,mommin,mommax,nbins,costmin,costmax,0.0,maxgap);
  TProfile2D* tsgapP = new TProfile2D("tsgapP","Theta scatter gap;Momentum (MeV/c);cos(#theta);gap (mm)",nbins,mommin,mommax,nbins,costmin,costmax,0.0,maxgap);
  TProfile2D* psgapP = new TProfile2D("psgapP","Phi scatter gap;Momentum (MeV/c);cos(#theta);gap (mm)",nbins,mommin,mommax,nbins,costmin,costmax,0.0,maxgap);
  TProfile2D* allgapP = new TProfile2D("allgapP","Total gap;Momentum (MeV/c);cos(#theta);gap (mm)",nbins,mommin,mommax,nbins,costmin,costmax,0.0,maxgap);
  mfgapP->SetStats(0);
  tsgapP->SetStats(0);
  psgapP->SetStats(0);
  allgapP->SetStats(0);
  mfgapP->SetMaximum(0.1);
  tsgapP->SetMaximum(2.0);
  psgapP->SetMaximum(2.0);
  allgapP->SetMaximum(2.0);
  mfgapP->SetMinimum(0.0);
  tsgapP->SetMinimum(0.0);
  psgapP->SetMinimum(0.0);
  allgapP->SetMinimum(0.0);
  TProfile2D* effP[4] = {mfgapP,tsgapP,psgapP,allgapP};

  // loop over samples
  for (unsigned isamp=0;isamp<nsamples;++isamp){
    double momval = gRandom->Uniform(mommin,mommax);
    // initial position is random inside range (+-)
    FourV pos;
    pos._x = gRandom->Uniform(-xrange,xrange);
    pos._y = gRandom->Uniform(-xrange,xrange);
    pos._z = gRandom->Uniform(-xrange,xrange);
    // initial time is arbitrary: set to 0
    pos._t = 0.0;
    FourV mom;
    // random momentum directions
    double cost = gRandom->Uniform(costmin,costmax);
    double sint = sqrt(1-cost*cost);
    double momphi = gRandom->Uniform(-pi,pi);
    mom._x = momval*sint*cos(momphi);
    mom._y = momval*sint*sin(momphi);
    mom._z = momval*cost;
    mom.SetEnergy(mass);
    // construct original helix from these
    HelixPars pars;
    HelixFromMom(pos,mom,charge,B,pars);
    // relativistic beta
    double beta = momval/mom._t;
// time range is assumed to be symmetric accross the measurement range.  Take
// the smaller of transverse or Z
    double tmax(10);
    double radius = pars_[rad_];
    double drange;
    if(2*radius < rmax) // spiral
      drange = zmax/cost;
    else
      drange = 2.0*radius*asin(rmax/(2.0*radius));
    double trange = drange/beta*c;
    trange = trange < tmax ? trange : tmax;
    // choose a random time for the effect
    FourV testpos;
    testpos._t = gRandom->Uniform(pars_[t0_]-trange,pars_[t0_]+trange);
    HelixPos(pars,testpos);
    // scattering sigma
    double scatsig = 14*sqrt(xfrac)/(momval*beta);
    // energy (momentum) loss fraction RMS
    double momfsig = 3.0*xfrac/(momval*beta*beta);
   // compute 1st order change in parameters
    HelixPars allpars = pars;
    // loop over scattering directions
    for(int ieff=0;ieff<3;++ieff){
    // choose correct sigma
      double sigma = ieff == 0 ? momfsig : scatsig;
      // compute random sample of this as a Gaussian (should be Moliere/landau!!, FIXME!!)
      double delta = gRandom->Gaus(0.0,sigma);
      HelixPars d1pars;
      HelixChange(pars,d1pars, testpos._t, delta, (effect)ieff);
      FourV d1pos;
      d1pos._t = testpos._t;
      HelixPos(d1pars,d1pos);
      // compute the gap between positions
      FourV gap = d1pos - testpos;
      effP[ieff]->Fill(momval,cost,gap.mag3());
      // accumulate the differences
      HelixPars apars = allpars;
      HelixChange(apars,allpars, testpos._t, delta, (effect)ieff);
    }
    // total
    FourV d1pos;
    d1pos._t = testpos._t;
    HelixPos(allpars,d1pos);
    // compute the gap between positions
    FourV allgap = d1pos - testpos;
     effP[3]->Fill(momval,cost,allgap.mag3());
  }
  // draw comparisons
  TCanvas* hacc = new TCanvas("hacc","Helix Accuracy",1000,600);
  hacc->Divide(2,2);
  hacc->cd(1);
  mfgapP->Draw("colorz");
  hacc->cd(2);
  tsgapP->Draw("colorz");
  hacc->cd(3);
  psgapP->Draw("colorz");
  hacc->cd(4);
  allgapP->Draw("colorz");
}

// change in residual for a small change in parameter.
void ParamDerivs(HelixPars const& pars, double time,double dvec[6]){
  // small increments for each parameter
  double dpar[6] = {0.1,0.1,0.1,0.1,0.001,0.1};
  // position at the test point
  FourV opos;
  opos._t = time;
  HelixPos(pars,opos);
// directions at test point
  FourV momdir;
  momdir._t = time;
  FourV thetadir;
  thetadir._t = time;
  HelixChangeDir(pars,theta1,thetadir);
  FourV phidir;
  phidir._t = time;
  HelixChangeDir(pars,theta2,phidir);
// now increment parameter
  HelixPars cpars = pars;
  for(unsigned ipar=0;ipar<5;++ipar){
     switch(ipar) {
      case 0:
	cpars_[rad_] += dpar[ipar];
	break;
      case 1:
	cpars_[lam_] += dpar[ipar];
	break;
      case 2:
	cpars_[cx_] += dpar[ipar];
	break;
      case 3:
	cpars_[cy_] += dpar[ipar];
	break;
      case 4:
	cpars_[phi0_] += dpar[ipar];
	break;
      case 5:
	cpars_[t0_] += dpar[ipar];
	break;
    }
    FourV cpos;
    cpos._t = time;
    HelixPos(cpars,cpos);
    FourV delta = opos-cpos;
// projection perpendicular to momentum
    double dphi = delta.dot3(phidir);
    double dtheta = delta.dot3(thetadir);
    // arbitrary sign
    double del = copysign(sqrt(dphi*dphi+dtheta*dtheta),dphi);
    // derivative
    dvec[ipar] = del/dpar[ipar];    
  }
  // t0 is a special case; derivative is constant, sign is arbitrary and unrelated to geometry (LR ambiguity).
  dvec[5] = copysign(vdrift,cos(time/1.0e-5));
}


void ParamCovar(double mommag, double mass, double cost, unsigned nmeasure, double trange, double xfrac, double msig) {
// create helix
  FourV pos(0.0,0.0,0.0,0.0);
  double sint = sqrt(1.0 - cost*cost);
  FourV mom(mommag*sint,0.0,mommag*cost,0.0);
  mom.SetEnergy(mass);
  HelixPars pars;
  HelixFromMom(pos, mom, -1.0, 1.0, pars);
  // scattering per measrurement
  double xfmeas = xfrac/nmeasure;
  // scattering sigma per measurement
  double beta = pars.beta();
  double scatsig = 14.0*sqrt(xfmeas)/(mommag*beta);
  // measurement weight
  double wm = 1.0/msig;
  // time step per measurement; time assumed symmetric
  double tstep = 2*trange/nmeasure;
  // initialize covariance
// loop over measurements
  for(unsigned imeas = 0;imeas<nmeasure;++imeas){
    double time = imeas*tstep-trange;
// compute weight for this measurement
  
    
// scattering derivatives
    HelixPars dt1pars;
    HelixChange(pars,dt1pars,time,scatsig,theta1);
    HelixPars dt2pars;
    HelixChange(pars,dt2pars,time,scatsig,theta2);
// update covariance


  }
}

void ParamCovarTest(double minmom, double maxmom, double mass, double mincost, double maxcost, unsigned nbins, unsigned nmeasure,double trange, double xfrac, double msig) {
  double momstep = (maxmom-minmom)/(nbins-1);
  double coststep = (maxcost-mincost)/(nbins-1);
  // loop over momenta
  for(unsigned mombin=0;mombin<nbins;++mombin){
    double momval = minmom + mombin*momstep;
    for(unsigned costbin=0;costbin<nbins;++costbin){
      double cost = mincost + costbin*coststep;
// compute covariance for these values
      ParamCovar(momval,mass,cost,nmeasure,trange,xfrac,msig);
    }
  }
}

void ComputeDOCA(DHit const& hit, HelixPars const& hpars, HDOCA& hdoca) {
  PVec const& pars = hpars_;
  // compute helix azimuth angle from the hit z position
  double phi = hit._wpos._z/pars[lam_] + pars[phi0_];
  // denomiator
  // angle difference
  double phibar = phi-hit._weta;
  double sinphibar = sin(phibar);
  double cosphibar = cos(phibar);
  double sineta = sin(hit._weta);
  double coseta = cos(hit._weta);
  double l2 = pars[lam_]*pars[lam_];
  double r2 = pars[rad_]*pars[rad_];
  double s2 = sinphibar*sinphibar;
  double m2 = hpars._mbar*hpars._mbar;
  double denom =  l2 + s2*r2;
  double Factor = pars[lam_]/sqrt(denom);
  double dx = pars[cx_] - hit._wpos._x;
  double dy = pars[cy_] - hit._wpos._y;
  double ddot = -sineta*dx + coseta*dy;
  // doca; this approximation ignores 2nd order terms.  It becomes
  // innacurate when tandip is small and z large
  double doca = -Factor*(pars[rad_]*cosphibar - ddot );
//  cout << " doca " << " Factor " << Factor << " ddot " << ddot << " r*cosphi " << pars[rad_]*cosphibar << endl;
  hdoca._doca = doca;
  // propagation time
  double tprop = hpars.time(hit._wpos._z);
  // time estimate for this hit.  This is equivalent to DOCA in the time domain
  hdoca._toca = tprop + fabs(doca)/vdrift;
//  cout << "wire t = " << hit._wpos._t << " hit t " << doca._toca <<  " tprop " << tprop << " htime " << hit._hpos._t <<  " hit z " << hit._wpos._z << endl;
  double ambig = doca > 0.0 ? 1.0 : -1.0;
  if(fabs(doca) > minDOCA)
    hdoca._ambig = ambig;
  else
    hdoca._ambig = 0.0;
  //  cout << "tresid = " << dt << " tpull = " << dt/hit._sigt << endl;
  // derivatives of DOCA
  PVec ddoca;
  ddoca[t0_] = 0.0; // no t0 dependence, purely geometric
  ddoca[cx_] = -Factor*sineta;
  ddoca[cy_] = Factor*coseta;
  // components that depend on phi; this ignores terms of magnitude DOCA/Radius
  ddoca[phi0_] = Factor*pars[rad_]*sinphibar;
  ddoca[rad_] = -Factor*(l2*cosphibar + pars[rad_]*s2*ddot)/denom;
    // this becomes inaccurate when z=0, but should have no effect on the fit
  ddoca[lam_] = -Factor*pars[rad_]*sinphibar*hit._wpos._z/l2; 
  // the extra terms are negligible, relative order DOCA/R and can be ignored
  //ddoca[lam_] = Factor*sinphibar*hit._wpos._z/l2 + sinphibar*(sinphibar + cosphibar*hit._wpos._z/(pars[rad_]*pars[lam_]))*doca/(pars[lam_]*denom);
  hdoca._ddoca = ddoca;
  // drift time residual.  This is signed by the ambiguity.  Use the true ambiguity
  // here, but in a real experiment it would need to be determined from data
  PVec dtdrift = ambig*(1.0/vdrift)*ddoca;
  // propagation time residual
  PVec dtprop;
  // add explicit dependence on parameters due to propagation
  dtprop[t0_] = 1.0; // direct dependence on t0
//  double f2 = -hit._wpos._z/(c*pars[lam_]*ebar);
  // the following terms are negligible and contribute and can be ignored
//  dtprop[lam_] = f2*(pars[rad_]*pars[rad_] + mbar*mbar)/(c*pars[lam_]);
//  dtprop[rad_] = f2*mbar*mbar/pars[rad_];
  hdoca._dtoca = dtprop + dtdrift;
//  cout << " dtprop[tlambda] " << dtprop[lam_] << " dtdrift[tlambda] " << dtdrift[lam_]
//  << " dtprop[rad] " << dtprop[rad_] << " dtdrift[rad] " << dtdrift[rad_] << endl;

}

void GenerateTimeHit(HelixPars const& pars, THit& thit, double tz) {
// compute the time for this hit location
  thit._hpos._t = pars.time(tz);
// find helix position and momentum at this z
  HelixPos(pars,thit._hpos);
// smear the time
  thit._t = gRandom->Gaus(thit._hpos._t,thit._sigt);
  // derivatives are trivial
  thit._dtoca[t0_] = 1.0;
}
 
void GenerateDriftHits(HelixPars const& pars, std::vector<DHit>& hits, unsigned nz, double deltaz, double sigdrift=0.2, double rmax=2.5) {
  hits.clear();
  double pbar = pars.pbar();
  double ebar = pars.ebar();
  double beta = pbar/ebar;
  double omega = pars.omega();
//  cout << "omega " << omega << "omega * lambda " << omega * pars_[lam_] << endl;
  // divide z range around 0
  for(int iz=-nz; iz<=(int)nz; ++iz){
    FourV hpos;
    double z = iz*deltaz;
    // compute the time from this position
    hpos._t = pars.time(z);
// find helix position and momentum at this z
    HelixPos(pars,hpos);
    double rho = sqrt(hpos._x*hpos._x+hpos._y*hpos._y);
//    if(rho > rmin && rho < rmax) {
      FourV mom;
      HelixMom(pars, hpos._t, mom);
      // generate a random azimuth angle for the wire
      double eta = gRandom->Uniform(-pi,pi);
      // define wire direction 
      FourV wdir(cos(eta),sin(eta),0.0,0.0);
      // define perp direction to wire and track
      FourV udir;
      CrossProduct(wdir,mom,udir);
      udir.normalize();
      // generate random drift; this is along the U direction
      double rdrift = gRandom->Uniform(-rmax,rmax); // sign defines ambiguity
      // define wire position
      FourV wpos = hpos;
      wpos._x -= rdrift*udir._x;
      wpos._y -= rdrift*udir._y;
      wpos._z -= rdrift*udir._z;
      double tdrift = fabs(rdrift/vdrift);
      wpos._t = hpos._t + tdrift;
      // create hit
      DHit hit;
      hit._hpos = hpos;
      hit._wpos = wpos;
      hit._weta = eta;
      hit._rdrift = rdrift;
      hit._sigt = sigdrift/vdrift;
      // generate random time smearing
      hit._t = gRandom->Gaus(wpos._t,hit._sigt);
      hits.push_back(hit);
//    }
  }
}


void PrintHit(DHit const& hit ) {
  cout << "Wire x " << hit._wpos._x << " Helix x " << hit._hpos._x 
    << " Wire y " << hit._wpos._y << " Helix y " << hit._hpos._y 
    << " Wire z " << hit._wpos._z << " Helix z " << hit._hpos._z 
    << " hit t " << hit._t << " helix t " << hit._hpos._t 
    << " eta " << hit._weta  << " rdrift " << hit._rdrift << endl;
}

void PrintDOCA(HDOCA const& doca ) {
  cout << "DOCA " << doca._doca << " TOCA " << doca._toca << endl;
}


void PrintHits(std::vector<DHit>const& hits) {
  for(auto hit: hits) {
    PrintHit(hit);
  }
}

void TestHelixFitDerivs(double cost=0.7, double momval=105.0,double charge=1.0,double mass=0.5,double deltaz=75.0) {
  double zrange = 3000.0;
  double momphi = pi/4.0;
  // parameter change for derivatives
  vector<double> delpars { 0.1, 0.025, 0.1, 0.1, 0.001, 0.1}; // small parameter changes for derivative calcs

  // position and momentum vectors
  FourV pos;
  FourV mom;
  double sint = sqrt(1-cost*cost);
  mom._x = momval*sint*cos(momphi);
  mom._y = momval*sint*sin(momphi);
  mom._z = momval*cost;
  mom.SetEnergy(mass);
// construct the helix from these
  HelixPars pars;
  HelixFromMom(pos,mom,charge,B,pars);
  pars_[t0_] = 2.8173; // random
  cout << "helix pars " << pars_ << " mass " << pars._mass << " mbar " << pars._mbar << endl;
// generate hits
  std::vector<DHit> hits;
  double pbar = pars.pbar();
  double ebar = pars.ebar();
  double beta = pbar/ebar;
  unsigned nz = ceil(zrange/(2*deltaz));
  GenerateDriftHits(pars,hits,nz,deltaz);
// print them
//  PrintHits(hits);
  // test derivatives; first setup parameter changes
  // positive change in helix
  vector<HelixPars> dpars(6,pars);
  for(size_t ipar=0;ipar < _npars;ipar++)
    dpars[ipar]._pars[ipar] += delpars[ipar];
// histogram differences
  vector<TH1F*> dderivh(_npars);
  vector<TGraph*> dderivg(_npars);
  vector<TH1F*> dtderivh(_npars);
  vector<TGraph*> dtderivg(_npars);
  string htitle, gtitle, hname, gname;
  for(int ipar=0;ipar<_npars;++ipar){
    hname = string("dderiv_")+string(parnames[ipar]) + string("_H");
    htitle = string("DOCA Change - Derivative*#Delta ")+string(parnames[ipar])+string(";#Delta DOCA (mm)");
    gname = string("dderiv_")+string(parnames[ipar]) + string("_G");
    gtitle = string("Derivative Change vs DOCA change: ") + string(parnames[ipar]) + string(";#Delta DOCA (mm);#Delta Dervi*#Delta")+string(parnames[ipar]) +string("(mm)");
    dderivh[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-0.01,0.01);
    dderivg[ipar] = new TGraph(hits.size());
    dderivg[ipar]->SetTitle(gtitle.c_str());
// time deriv
    hname = string("dtderiv_")+string(parnames[ipar]) + string("_H");
    htitle = string("Time Change - Time Derivative*#Delta ")+string(parnames[ipar])+string(";#Delta Time (ns)");
    gname = string("dtderiv_")+string(parnames[ipar]) + string("_G");
    gtitle = string("Time Derivative Change vs Time change: ") + string(parnames[ipar]) + string(";#Delta Time (ns);#Delta TDervi*#Delta")+string(parnames[ipar]) +string("(ns)");
    dtderivh[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-0.1,0.1);
    dtderivg[ipar] = new TGraph(hits.size());
    dtderivg[ipar]->SetTitle(gtitle.c_str());
  }
// time residuals, etc
  TH1F* tresid = new TH1F("tresid","Time Residuals;#Delta t(ns)",100,-25.0,25.0);
  TH1F* tpull = new TH1F("tpull","Time Pulls",100,-10.0,10.0);
  TH1F* tdrift = new TH1F("tdrift","Drift Time;ns",100,-1.0,50.0);
  TH1F* tres = new TH1F("tres","Time Resolution;ns",100,-0.2,0.2);
  unsigned ihit(0);
  for(auto hit : hits) {
    HDOCA hdoca; // nominal doca
    ComputeDOCA(hit,pars,hdoca);
//    cout << "doca " << hdoca._doca << " rdrift " << hit._rdrift << endl;
    tresid->Fill(hdoca._toca-hit._t);
    tpull->Fill((hdoca._toca-hit._t)/hit._sigt);
    tdrift->Fill(hit._wpos._t - hit._hpos._t);
    tres->Fill(hdoca._toca-hit._wpos._t);
    // compute DOCA for original and modified paeters 
    for(int ipar=0;ipar<_npars;++ipar){
      HDOCA hddoca;
      ComputeDOCA(hit,  dpars[ipar], hddoca);
      double ddoca = (hddoca._doca - hdoca._doca);
      double dderiv = hdoca._ddoca[ipar]*delpars[ipar];
//      cout << parnames[ipar] << " Hit z " << hit._wpos._z << " DeltaDoca " << ddoca << " dderiv " << dderiv << endl;
      dderivh[ipar]->Fill(ddoca - dderiv);
      dderivg[ipar]->SetPoint(ihit++,ddoca,dderiv);
// now test time derivatives
      double dtime = (hddoca._toca - hdoca._toca);
      double dtderiv = hdoca._dtoca[ipar]*delpars[ipar];
//      cout << parnames[ipar] << " Hit z " << hit._wpos._z << " D Hit Time " 
//	<< dtime << " dtderiv " << dtderiv 
//	<< endl;
      dtderivh[ipar]->Fill(dtime - dtderiv);
      dtderivg[ipar]->SetPoint(ihit++,dtime,dtderiv);
//      if(fabs(dtime-dtderiv) > 0.1) {
//	cout << "outlier par " << ipar << " dtime " << dtime << " dtderiv " << dtderiv 
//	<< " derivatives " << doca._dtoca << endl;
//	cout <<" original hit " << endl;
//	PrintHit(hit);
//	cout <<" modified hit " << endl;
//	PrintHit(dhit);
//      }
     }
  }

  TCanvas* dderivhc = new TCanvas("dderivhc","dderivhc",800,600);
  dderivhc->Divide(3,2);
  for(int ipar=0;ipar<_npars;++ipar){
    dderivhc->cd(ipar+1);
    dderivh[ipar]->Draw();
  };
  TCanvas* dderivgc = new TCanvas("dderivgc","dderivgc",800,600);
  dderivgc->Divide(3,2);
  for(int ipar=0;ipar<_npars;++ipar){
    dderivgc->cd(ipar+1);
    dderivg[ipar]->Draw("A*");
  };

  TCanvas* dtderivhc = new TCanvas("dtderivhc","dtderivhc",800,600);
  dtderivhc->Divide(3,2);
  for(int ipar=0;ipar<_npars;++ipar){
    dtderivhc->cd(ipar+1);
    dtderivh[ipar]->Draw();
  };
  TCanvas* dtderivgc = new TCanvas("dtderivgc","dtderivgc",800,600);
  dtderivgc->Divide(3,2);
  for(int ipar=0;ipar<_npars;++ipar){
    dtderivgc->cd(ipar+1);
    dtderivg[ipar]->Draw("A*");
  };

  TCanvas* residcan = new TCanvas("residcan","residcan",800,800);
  residcan->Divide(2,2);
  residcan->cd(1);
  tresid->Draw();
  residcan->cd(2);
  tpull->Draw();
  residcan->cd(3);
  tdrift->Draw();
  residcan->cd(4);
  tres->Draw();
}

void FitHelix(std::vector<DHit>const& hits, THit const& thit, HelixPars const& refpars, PVec& fpar, PMat& fcov, double& chisq) {

  // deweight input covarance
  fcov *= 1000000.0;
  //  cout <<" initial covariance " << endl << cov << endl;
  // invert this for the initial weight
  PMat fweight = fcov;
  bool inverse = fweight.Invert();
  if(!inverse){
    cout << "inversion of initial covariance failed " << endl;
    cout << fcov << endl;
    return;
  }
  //  cout <<" initial weight " << endl << fweight << endl;
  // input parameters are reference
  // initial weight vector
  PVec fbeta = fweight*fpar;
  //  cout << " initial beta " << endl << fbeta << endl;
  // now process the hits
  SMatrix<double, 1,1, MatRepSym<double,1> > I1 = SMatrixIdentity();
  PDer pder;
  chisq = 0.0;
  for(auto hit : hits) {
    // compute residual WRT the reference parameters
    HDOCA refdoca;
    ComputeDOCA(hit,refpars, refdoca);
    // convert the parameter derivatives and residual into a weight and weight constraint
    PVec hbeta;
    PMat hweight;
    double hwt;
    double resid;
    PVec pderivs;
    if(refdoca._ambig != 0){
      // ambiguity is defined so use time residual (drift)
      hwt = 1.0/hit._sigt;
      resid = (hit._t - refdoca._toca);
      pderivs = refdoca._dtoca;
    } else {
      // too close to the wire to use time residual: use DOCA to the wire as the measurement
      // error is approximate in this case
      // double because of sign of DOCA
      static double sqrt12 = sqrt(12.0);
      hwt = sqrt12/(2.0*minDOCA);
      resid = -refdoca._doca; // 'measured' value is 0
      pderivs = refdoca._ddoca;
    }
    // weight matrix and weight vector for this hit
    // must convert manually from SVector SMatrix[N][1], ugh!
    for(size_t ipar=0; ipar<_npars; ++ipar) {
      pder[ipar][0] = pderivs[ipar]*hwt;
    }
    //      cout << doca._dtoca << endl <<  pder << endl;
    hweight = Similarity(pder,I1);
    //	cout << "hweight " << endl << hweight << endl;
    // reference beta for this hit
    PVec hbetaref = hweight*fpar;
      // change in parameters from this residual
    PVec delta = pderivs*resid*hwt*hwt;
      //	PVec delta = gendoca._dtoca*tresid*hwt*hwt;
      //      cout << "hweight " << hweight << endl;
      //      cout << "ref =   " << hbetaref << endl;
      //      cout << "delta   " << delta << endl;
      // add change WRT reference; sign convention reflects resid = measure - prediction
    hbeta = hbetaref + delta; // check sign
    // compute chisquared contribution for this hit
    // use covariance from previous hit
    PMat prevcov = fweight;
    if(prevcov.Invert()) {
    // compute chisquared increment for this hit
      PVec prevpar = prevcov*fbeta;
      double parerr2 =  Similarity(prevcov,pderivs);
      double hiterr2 = 1.0/(hwt*hwt);
      double chierr2 = parerr2 + hiterr2;
      double dresid =  Dot(pderivs,prevpar-fpar); 
      double chival = resid - dresid;
      double dchisq = chival*chival/chierr2;
//      cout << " resid " << resid << " dresid " << dresid
//      << " chival " << chival << " parerr " << sqrt(parerr2)
//      << " hiterr " << sqrt(hiterr2) << " chierr2 " << sqrt(chierr2) 
//      << " dchisq " << dchisq << " chisq " << chisq <<  endl;
      chisq += dchisq;
    }
    // update the fit for this hit
    fbeta += hbeta;
    fweight += hweight;
  }
  // add time hit effects
  double hwt = 1.0/thit._sigt;
  double resid = thit._t - refpars.time(thit._hpos._z);
  for(size_t ipar=0; ipar<_npars; ++ipar) {
    pder[ipar][0] = thit._dtoca[ipar]*hwt;
  }
  PMat hweight = Similarity(pder,I1);
  PVec hbetaref = hweight*fpar;
  PVec delta = thit._dtoca*resid*hwt*hwt;
  PVec hbeta = hbetaref + delta; // check sign
  PMat prevcov = fweight;
  if(prevcov.Invert()) {
    // compute chisquared increment for this hit
      PVec prevpar = prevcov*fbeta;
      double parerr2 =  Similarity(prevcov,thit._dtoca);
      double hiterr2 = 1.0/(hwt*hwt);
      double chierr2 = parerr2 + hiterr2;
      double dresid =  Dot(thit._dtoca,prevpar-fpar); 
      double chival = resid - dresid;
      double dchisq = chival*chival/chierr2;
//      cout << " resid " << resid << " dresid " << dresid
//      << " chival " << chival << " parerr " << sqrt(parerr2)
//      << " hiterr " << sqrt(hiterr2) << " chierr2 " << sqrt(chierr2) 
//      << " dchisq " << dchisq << " chisq " << chisq <<  endl;
      chisq += dchisq;
  }
  fbeta += hbeta;
  fweight += hweight;
  // compute chisquared contribution for this hit
  // use covariance from previous hit

  //  cout << " final beta " << endl << fbeta << endl;
  //  cout << " final weight " << endl <<  fweight << endl;
  // recover the fit parameters
  fcov = fweight;
  inverse = fcov.Invert();
  if(!inverse){
    cout << "inversion of final weight failed " << endl;
    cout << fweight << endl;
    return;
  }
  // update parameters for output
  fpar = fcov*fbeta;
}

void TestHelixFit(unsigned ntries, unsigned maxniter=10, double reffac=3.0, double momval=105.0,double charge=1.0,double genmass=0.5,double fitmass=0.5,double deltaz=75.0,double sigdrift=0.2, double sigtime=0.3, bool invert=false) {
  double zrange = 3000.0;
  unsigned nz = ceil(zrange/(2*deltaz));
  cout << "will generate " << 2*nz+1 << " z planes " << endl;
  double scale = 1.0/sqrt(nz);
  double costmin(0.5);
  double costmax(0.7);
  vector<double> parsig(6);
  parsig[rad_] = 0.5*scale; // approximate error in radius 
  parsig[cx_]= 0.5*scale; // center x
  parsig[cy_]= 0.5*scale; // in center y
  parsig[lam_] = 0.1*scale; // in tan lambda
  parsig[phi0_]= 0.005*scale; // in phi0
  parsig[t0_]= 5.0*scale; // in t0
  vector<TGraph*> dpref(_npars);
  vector<TGraph*> dpgen(_npars);
  vector<TH1F*> dprefh(_npars);
  vector<TH1F*> dpgenh(_npars);
  vector<TH1F*> dpullrefh(_npars);
  vector<TH1F*> dpullgenh(_npars);
  vector<TH1F*> fiterrh(_npars);
  string gtitle, htitle, hname;
  for(int ipar=0;ipar<_npars;++ipar){
//    gtitle = string("Fit vs Ref: ") + string(parnames[ipar]) + string(";Ref ;Fit");
//    dpref[ipar] = new TGraph(ntries);
//    dpref[ipar]->SetTitle(gtitle.c_str());
//    dpref[ipar]->SetMarkerStyle(20);
    gtitle = string("Fit vs Gen: ") + string(parnames[ipar]) + string(";Gen")
    + parunits[ipar] + string(";Fit") + parunits[ipar];
    dpgen[ipar] = new TGraph(ntries);
    dpgen[ipar]->SetTitle(gtitle.c_str());
    dpgen[ipar]->SetMarkerStyle(20);
//    htitle = string("#Delta")+string(parnames[ipar]) + string("Fit - Ref");
//    hname = string("D")+string(parnames[ipar])+string("FR");
//    dprefh[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-10*parsig[ipar],10*parsig[ipar]);
    htitle = string("#Delta")+parnames[ipar] + string("Fit - Gen")+string(";") + parnames[ipar] + parunits[ipar];
    hname = string("D")+parnames[ipar]+string("FG");
    dpgenh[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-5*parsig[ipar],5*parsig[ipar]);
//    htitle = string("#Delta")+string(parnames[ipar]) + string("Pull Fit - Ref");
//    hname = string("DP")+string(parnames[ipar])+string("FR");
//    dpullrefh[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-10,10);
    htitle = string("#Delta")+parnames[ipar] + string("Pull Fit - Gen");
    hname = string("DP")+parnames[ipar]+string("FG");
    dpullgenh[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,-10,10);
    htitle = parnames[ipar] + string("Fit Error;#sigma") + parnames[ipar] + parunits[ipar];
    hname = parnames[ipar]+string("FE");
    fiterrh[ipar] = new TH1F(hname.c_str(),htitle.c_str(),100,0.0,2*parsig[ipar]);
  }
  TH1F* chisqh = new TH1F("chisqh","Chisquared",200,0.0,4.0*nz);
  TH1F* chisqndofh = new TH1F("chisqndofh","Chisquared/DOF",200,0.0,5.0);
  TH1F* chisqprobh = new TH1F("chisqprobh","Chisquared Probability",200,0.0,1.0);
  TH1F* logchisqprobh = new TH1F("logchisqprobh","log10 (Chisquared Probability)",200,-50.0,0.0);
  // correlation matrix histogram
  TH2F* corravg = new TH2F("corravg","Average correlation matrix magnitudes",_npars,-0.5,_npars-0.5,_npars, -0.5,_npars-0.5);
// momentum plots
  TH1F* mommagh = new TH1F("mommagh","Fit Momentum",100,momval-0.5,momval+0.5);
  TH1F* momresh = new TH1F("momresh","Fit Momentum Resolution",100,-0.5,0.5);
  TH1F* momerrh = new TH1F("momerrh","Fit Momentum Error",100,0,0.5);
  TH1F* mompullh = new TH1F("mompullh","Fit Momentum Pull",100,-10.0,10.0);
  // loop over tries
  for(unsigned itry = 0; itry < ntries; itry++){
    double cost = gRandom->Uniform(costmin,costmax);
    double momphi = gRandom->Uniform(-pi,pi);
    //    cout << "momphi " << momphi << endl;
    // position and momentum vectors
    FourV pos;
    pos._x = gRandom->Uniform(-5,5);
    pos._y = gRandom->Uniform(-5,5);
    pos._z = gRandom->Uniform(-200,200);
    FourV mom;
    double sint = sqrt(1-cost*cost);
    mom._x = momval*sint*cos(momphi);
    mom._y = momval*sint*sin(momphi);
    mom._z = momval*cost;
    mom.SetEnergy(genmass);
    // construct the helix from these
    HelixPars genpars;
    HelixFromMom(pos,mom,charge,B,genpars);
    genpars_[t0_] = gRandom->Uniform(-10.0,10.0);
    // generate hits
    std::vector<DHit> hits;
    GenerateDriftHits(genpars,hits,nz,deltaz,sigdrift);
//    cout << "nz " << nz << " nhits " << hits.size() << endl;
//    cout << "generator pars " << endl << genpars_ << endl;
    THit thit;
    thit._sigt = sigtime;
    GenerateTimeHit(genpars,thit,timez);
    // construct time domain fit chisquared
    // initialize the covariance matrix from initial errors
    // randomize initial parameters
    PMat cov;
    HelixPars refpars(genpars);
    // update to the fit mass
    refpars._mass = fitmass;
    refpars._mbar = genpars._mbar*fitmass/genmass;
    for(unsigned ipar=0;ipar < _npars; ipar++){
      cov[ipar][ipar] = parsig[ipar];
      refpars_[ipar] += gRandom->Gaus(0.0,parsig[ipar]*reffac);
    }
    if(invert)refpars.invert();
    PVec fpar = refpars_;
    PMat fcov = cov;
//    cout << "ref pars " << iter << endl << fpar << endl;
    // iterate fit
    bool test(false);
    unsigned iter(0);
    double chisq(0.0);
    int ndof = hits.size()-6 + 1;
    while ((!test) && iter < maxniter) {
      FitHelix(hits,thit,refpars,fpar,fcov,chisq);
  //    cout << "pars after fit " << iter << endl << fpar << endl;
      refpars_ = fpar;
      iter++;
    }
//      cout << "N hits = " << hits.size() << endl;
    //  cout << "true parameters          " << pars_ << endl;
    //  cout << "fit reference parameters " << refpars_ << endl;
    //  cout << "fit result parameters    " << fpar << endl;
    //   cout << "fit chisq " << chisq << " ndof " << ndof << endl;
    PVec cdiag;
    for(size_t ipar=0;ipar<_npars;ipar++){
      cdiag[ipar] = sqrt(fcov[ipar][ipar]);
    }
    //  cout << "fit covar diag sqrt      " << cdiag << endl;
    //cout << "fit result covariance " << endl << fcov << endl;
    for(size_t ipar=0;ipar< _npars; ipar++){
//      dpref[ipar]->SetPoint(itry, refpars_[ipar], fpar[ipar]);
      dpgen[ipar]->SetPoint(itry, genpars_[ipar], fpar[ipar]);
//      dprefh[ipar]->Fill(fpar[ipar]-refpars_[ipar]);
      dpgenh[ipar]->Fill(fpar[ipar]-genpars_[ipar]);
//      dpullrefh[ipar]->Fill((fpar[ipar]-refpars_[ipar])/cdiag[ipar]);
      dpullgenh[ipar]->Fill((fpar[ipar]-genpars_[ipar])/cdiag[ipar]);
      fiterrh[ipar]->Fill(cdiag[ipar]);
    }
// accumulate average correlation matrix
    PMat corr = fcov;
    for(unsigned ipar=0; ipar <_npars;ipar++){
      for(unsigned jpar=ipar;jpar < _npars; jpar++){
	corr[ipar][jpar] /= cdiag[ipar]*cdiag[jpar];
	corravg->Fill(ipar,jpar,fabs(corr[ipar][jpar]));
      }
    }
//    cout << "Correlation matrix " << endl << corr << endl;
// accumulate chisquared info
    double chisqprob = TMath::Prob(chisq,ndof);
    chisqh->Fill(chisq);
    chisqndofh->Fill(chisq/ndof);
    chisqprobh->Fill(chisqprob);
    logchisqprobh->Fill(log10(chisqprob));

    // compute momentum and momentum error
    PVec momproj;
    double factor = refpars._mass/refpars._mbar;
    momproj[rad_] = factor*refpars_[rad_]/refpars.pbar();
    momproj[lam_] = factor*refpars_[lam_]/refpars.pbar();
    double momerr2 = Similarity(fcov,momproj);
//    double mommag = fabs(factor)*sqrt(fpar[rad_]*fpar[rad_] + fpar[lam_]*fpar[lam_] );
    FourV fmom;
    HelixMom(refpars,fpar[t0_],fmom);
    double mommag = fmom.mag3();
//cout << "mommag " << mommag << " vecmag " << fmom.mag3() << endl;
    // test
    mommagh->Fill(mommag);
    momresh->Fill(mommag - momval);
    momerrh->Fill(sqrt(momerr2));
    mompullh->Fill( (mommag-momval)/sqrt(momerr2) );


  }


//  TCanvas* dprefcan = new TCanvas("dprefcan","dprefcan",800,600);
//  dprefcan->Divide(3,2);
//  for(int ipar=0;ipar<_npars;++ipar){
//    dprefcan->cd(ipar+1);
//    dpref[ipar]->Draw("AP");
//  }
  TCanvas* dpgencan = new TCanvas("dpgencan","dpgencan",800,600);
  dpgencan->Divide(3,2);
  for(int ipar=0;ipar<_npars;++ipar){
    dpgencan->cd(ipar+1);
    dpgen[ipar]->Draw("AP");
  }
//  TCanvas* dprefhcan = new TCanvas("dprefhcan","dprefhcan",800,600);
//  dprefhcan->Divide(3,2);
//  for(int ipar=0;ipar<_npars;++ipar){
//    dprefhcan->cd(ipar+1);
//    dprefh[ipar]->Draw();
//  };
  TCanvas* dpgenhcan = new TCanvas("dpgenhcan","dpgenhcan",800,600);
  dpgenhcan->Divide(3,2);
  for(int ipar=0;ipar<_npars;++ipar){
    dpgenhcan->cd(ipar+1);
    dpgenh[ipar]->Fit("gaus");
  }
//  TCanvas* dpullrefhcan = new TCanvas("dpullrefhcan","dpullrefhcan",800,600);
//  dpullrefhcan->Divide(3,2);
//  for(int ipar=0;ipar<_npars;++ipar){
//    dpullrefhcan->cd(ipar+1);
//    dpullrefh[ipar]->Draw();
//  };
  TCanvas* dpullgenhcan = new TCanvas("dpullgenhcan","dpullgenhcan",800,600);
  dpullgenhcan->Divide(3,2);
  for(int ipar=0;ipar<_npars;++ipar){
    dpullgenhcan->cd(ipar+1);
    dpullgenh[ipar]->Fit("gaus");
  }
  TCanvas* fiterrhcan = new TCanvas("fiterrhcan","fiterrhcan",800,600);
  fiterrhcan->Divide(3,2);
  for(int ipar=0;ipar<_npars;++ipar){
    fiterrhcan->cd(ipar+1);
    fiterrh[ipar]->Draw();
  }
  TCanvas* corrcan = new TCanvas("corrcan","corrcan",600,600);
  corrcan->Divide(1,1);
  corrcan->cd(1);
  corravg->Scale(1.0/float(ntries));
  corravg->SetStats(0);
  gPad->SetLogz(); 
  TAxis* xax = corravg->GetXaxis();
  TAxis* yax = corravg->GetYaxis();
  for (unsigned ipar=0; ipar < _npars; ipar++){
    xax->SetBinLabel(ipar+1,parnames[ipar].c_str());
    yax->SetBinLabel(ipar+1,parnames[ipar].c_str());
  }
  corravg->Draw("colorztext0");
  // chisq
  TCanvas* chisqcan = new TCanvas("chisqcan","chisqcan",800,800);
  chisqcan->Divide(2,2);
  chisqcan->cd(1);
  chisqndofh->Draw();
  chisqcan->cd(2);
  chisqh->Draw();
  chisqcan->cd(3);
  chisqprobh->Draw();
  chisqcan->cd(4);
  logchisqprobh->Draw();
// momentum
  TCanvas* momcan = new TCanvas("momcan","momcan",800,800);
  momcan->Divide(2,2);
  momcan->cd(1);
  momresh->Fit("gaus");
  momcan->cd(2);
  mommagh->Draw();
  momcan->cd(3);
  momerrh->Draw();
  momcan->cd(4);
  mompullh->Fit("gaus");
}


