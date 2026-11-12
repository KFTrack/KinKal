//
// Test helix intersections with a surface the helix crosses, overshoots slightly and crosses back (KFTrack/KinKal#233).
// The query starts on the outgoing leg (outside the tolerance slab, or at the turning point) and ends after the
// return crossing, so intersect must return the return crossing, which is the only root inside the query.
// The first sign-change step straddles the turning point; the refinement must stay inside that step.
// Surfaces: a plane (Rectangle) and a cylinder parallel to the field; LoopHelix and CentralHelix, both time directions.
// Second part: a helix whose axis is tilted from the plane normal starts just in front of the plane, so the point on
// the helix axis at the query start can already be past the plane; the crossing just after the start must be found.
//
#include "KinKal/General/ParticleState.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/CentralHelix.hh"
#include "KinKal/Geometry/Rectangle.hh"
#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Geometry/Intersection.hh"
#include "KinKal/Geometry/Intersect.hh"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <memory>
#include <string>

using namespace KinKal;

enum class SurfType { plane, cylinder };

// startfrac: query start in half-widths before the turning point (0 = at it; > 1 = before the first crossing)
// expect the first crossing when the query starts before it, otherwise the return crossing
template <class HELIX> int surfaceReturn(SurfType stype, TimeDir tdir, double excursion, double tol, double startfrac, std::string const& name) {
  double sign = tdir == TimeDir::forwards ? 1.0 : -1.0;
  // 100 MeV/c transverse, 20 MeV/c along a 1 T field, muon mass; R = 333.6 mm
  ParticleState state(VEC3(0.0,0.0,0.0),VEC3(100.0,0.0,20.0),0.0,105.658,1);
  HELIX helix(state,VEC3(0.0,0.0,1.0),TimeRange(-100.0,100.0));
  double rbend = helix.bendRadius();
  double omega = helix.transverseSpeed()/rbend; // |dphi/dt|
  double tmax = sign*0.5*M_PI/omega; // the turning point: a quarter turn from the start
  std::unique_ptr<Surface> surf;
  double halfwidth(0.0); // time from a root to the turning point
  std::string sname;
  if(stype == SurfType::plane){
    // plane normal to the helix radius at tmax, excursion inside the helix circle
    double radius = std::fabs(helix.position3(tmax).X());
    halfwidth = std::acos(1.0-excursion/radius)/omega;
    surf = std::make_unique<Rectangle>(VEC3(1.0,0.0,0.0),VEC3(0.0,1.0,0.0),helix.position3(tmax-sign*halfwidth),2000.0,2000.0);
    sname = "Plane";
  } else {
    // cylinder along the field whose axis is offset from the helix axis by ecc towards the helix position at tmax:
    // the distance from the cylinder axis is smallest (rbend-ecc) at tmax, and the helix dips excursion inside it
    double ecc = 0.5*rbend; // the same half-width as the plane to second order
    auto rdir = helix.position3(tmax) - helix.center(tmax);
    rdir.SetZ(0.0);
    auto cpos = helix.center(tmax) + ecc*rdir.Unit();
    cpos.SetZ(0.0);
    double crad = rbend - ecc + excursion;
    halfwidth = std::acos((rbend*rbend + ecc*ecc - crad*crad)/(2.0*rbend*ecc))/omega;
    surf = std::make_unique<Cylinder>(VEC3(0.0,0.0,1.0),cpos,crad,2000.0);
    sname = "Cylinder";
  }
  double first = tmax - sign*halfwidth; // the crossing into the excursion
  double ret = tmax + sign*halfwidth; // the return crossing
  double expected = startfrac > 1.0 ? first : ret;
  double qstart = tmax - sign*startfrac*halfwidth;
  double qend = tmax + sign*1.5*halfwidth;
  TimeRange query = tdir == TimeDir::forwards ? TimeRange(qstart,qend) : TimeRange(qend,qstart);
  auto inter = KinKal::intersect(helix,*surf,query,tol,tdir);
  bool ok = inter.good() && std::fabs(inter.time_-expected)*helix.speed(expected) < tol;
  std::cout << std::setprecision(12) << name << " " << sname << (tdir == TimeDir::forwards ? " forwards" : " backwards")
    << " start " << startfrac << " half-widths before the turn, excursion " << excursion << " mm tol " << tol << " mm\n"
    << "  query [" << query.begin() << ", " << query.end() << "] ns, start "
    << surf->distance(helix.position3(qstart)) << " mm from the surface\n"
    << "  roots: first " << first << ", return " << ret << ", expected " << expected << "\n"
    << "  intersect: time " << inter.time_ << " onsurface " << inter.onsurface_ << " inbounds "
    << inter.inbounds_ << " inrange " << inter.inrange_ << "  -> " << (ok ? "OK" : "FAIL") << "\n";
  return ok ? 0 : 1;
}

// a 35.5 MeV/c electron in a 2 T field tilted by alpha from the plane normal (z) starts dist0 in front of the plane,
// moving towards it with cos = -0.53 to the normal; the helix loops within the query; beta is the momentum azimuth
template <class HELIX> int axisTilt(TimeDir tdir, double alpha, double beta, std::string const& name) {
  double sign = tdir == TimeDir::forwards ? 1.0 : -1.0;
  double tol(0.1), dist0(0.48), cost(-0.53), pmom(35.5), qlen(0.25);
  double sint = std::sqrt(1.0-cost*cost);
  // in the search direction the particle moves towards the plane: reverse the momentum when searching backwards
  VEC3 mom = sign*pmom*VEC3(sint*std::cos(beta),sint*std::sin(beta),cost);
  ParticleState state(VEC3(0.0,0.0,0.0),mom,0.0,0.511,-1);
  VEC3 bnom(2.0*std::sin(alpha),0.0,2.0*std::cos(alpha));
  HELIX helix(state,bnom,TimeRange(-100.0,100.0));
  Rectangle plane(VEC3(0.0,0.0,1.0),VEC3(1.0,0.0,0.0),VEC3(0.0,0.0,-dist0),2000.0,2000.0);
  TimeRange query = tdir == TimeDir::forwards ? TimeRange(0.0,qlen) : TimeRange(-qlen,0.0);
  // the expected crossing, by bisection over the query (the plane distance is monotonic over it)
  double ta(0.0), tb(sign*qlen);
  double fa = plane.distance(helix.position3(ta));
  if(fa*plane.distance(helix.position3(tb)) >= 0.0){
    std::cout << name << " axis tilt: no crossing in the query, bad test geometry" << std::endl;
    return 1;
  }
  for(unsigned i=0; i < 100; ++i){
    double tm = 0.5*(ta+tb);
    if(plane.distance(helix.position3(tm))*fa > 0.0) ta = tm; else tb = tm;
  }
  double expected = 0.5*(ta+tb);
  auto inter = KinKal::intersect(helix,plane,query,tol,tdir);
  bool ok = inter.good() && std::fabs(inter.time_-expected)*helix.speed(expected) < tol;
  auto ray = helix.linearize(0.0);
  std::cout << std::setprecision(12) << name << " axis tilt " << alpha << " azimuth " << beta
    << (tdir == TimeDir::forwards ? " forwards" : " backwards") << " query [" << query.begin() << ", " << query.end()
    << "] ns, radius " << helix.bendRadius() << " mm, linearized ray start " << plane.distance(ray.start())
    << " mm from the plane\n  expected " << expected << ", intersect: time " << inter.time_ << " onsurface "
    << inter.onsurface_ << " inbounds " << inter.inbounds_ << " inrange " << inter.inrange_ << "  -> " << (ok ? "OK" : "FAIL") << "\n";
  return ok ? 0 : 1;
}

int main() {
  int nfail(0), ncase(0);
  for(auto stype : {SurfType::plane, SurfType::cylinder}){
    for(auto tdir : {TimeDir::forwards, TimeDir::backwards}){
      // tol = 0.1 mm: the first sign-change step straddles the turning point
      // startfrac 0.5: the query starts outside the slab on the outgoing leg; 0.0: at the turning point
      for(double startfrac : {0.5, 0.0}){
        nfail += surfaceReturn<LoopHelix>(stype,tdir,0.15,0.1,startfrac,"LoopHelix"); ++ncase;
        nfail += surfaceReturn<CentralHelix>(stype,tdir,0.15,0.1,startfrac,"CentralHelix"); ++ncase;
      }
      // controls: tol = 0.05 mm (the sign-change step starts past the turn), and a query that contains both roots
      nfail += surfaceReturn<LoopHelix>(stype,tdir,0.15,0.05,0.5,"control LoopHelix"); ++ncase;
      nfail += surfaceReturn<CentralHelix>(stype,tdir,0.15,0.05,0.5,"control CentralHelix"); ++ncase;
      nfail += surfaceReturn<LoopHelix>(stype,tdir,0.15,0.1,1.5,"control LoopHelix"); ++ncase;
      nfail += surfaceReturn<CentralHelix>(stype,tdir,0.15,0.1,1.5,"control CentralHelix"); ++ncase;
    }
  }
  int nfailaxis(0), ncaseaxis(0);
  for(auto tdir : {TimeDir::forwards, TimeDir::backwards}){
    for(double alpha : {0.0, 0.05}){ // alpha = 0: the axis point and the particle are at the same distance from the plane
      for(unsigned ibeta=0; ibeta < 8; ++ibeta){
        double beta = ibeta*M_PI/4.0;
        nfailaxis += axisTilt<LoopHelix>(tdir,alpha,beta,"LoopHelix"); ++ncaseaxis;
        nfailaxis += axisTilt<CentralHelix>(tdir,alpha,beta,"CentralHelix"); ++ncaseaxis;
      }
    }
  }
  std::cout << nfail << " of " << ncase << " return-crossing cases fail; " << nfailaxis << " of " << ncaseaxis << " axis-tilt cases fail" << std::endl;
  nfail += nfailaxis;
  return nfail == 0 ? 0 : 1;
}
