//
// Test the projected spatial covariance of a trajectory onto the local plane at an intersection (Intersection::pcov_)
// by comparing it with a toy MC: the trajectory parameters are randomized according to their covariance, the
// intersection is recomputed, and the spread of the intersection points in the local plane UV system is measured.
// Also tests consistency of pcov_ with KTRAJ::planeCovariance, that the covariance is filled for
// intersections that are on the surface but out of bounds, and that it is only filled on request.
// original author: David Brown (LBNL) 2026
//
#include "KinKal/Geometry/Intersection.hh"
#include "KinKal/Geometry/Intersect.hh"
#include "KinKal/Geometry/Plane.hh"
#include "KinKal/Geometry/Rectangle.hh"
#include "KinKal/Geometry/Cylinder.hh"
#include <iostream>
#include <string>
#include <random>
#include <cmath>
#include <getopt.h>

using namespace KinKal;
using std::cout;
using std::endl;

static struct option long_options[] = {
  {"ntoy",       required_argument, 0, 'n'  },
  {"seed",       required_argument, 0, 's'  },
  {"verbose",    required_argument, 0, 'v'  },
  {NULL, 0,0,0}
};

void print_usage() {
  printf("Usage: IntersectionCovariance --ntoy (number of toy MC trials) --seed (random seed) --verbose (verbosity)\n");
}

// set a diagonal parameter covariance
template <class KTRAJ> void setCovariance(KTRAJ& ktraj, std::array<double,6> const& sigmas) {
  for(size_t ipar=0;ipar < NParams(); ++ipar){
    for(size_t jpar=0;jpar < NParams(); ++jpar){
      ktraj.params().covariance()(ipar,jpar) = ipar == jpar ? sigmas[ipar]*sigmas[ipar] : 0.0;
    }
  }
}

// INTERF is a callable (KTRAJ const&, SURF const&, TimeRange const&, bool pcov) -> Intersection
template <class KTRAJ, class SURF, class INTERF> int testCovariance(std::string const& name, KTRAJ const& ktraj, SURF const& surf,
    TimeRange const& trange, INTERF interf, bool expectinbounds, unsigned ntoy, unsigned seed, int verbose) {
  static const double vartol(0.05); // relative tolerance on toy variance (statistical precision ~1% for 20K toys)
  static const double rhotol(0.03); // absolute tolerance on the toy correlation coefficient
  static const double ctol(1.0e-9); // relative tolerance on consistency checks
  auto inter = interf(ktraj,surf,trange,true);
  if(verbose > 0) cout << name << " " << inter << endl;
  if(!inter.onsurface_){
    cout << name << ": intersection not on surface" << endl;
    return -1;
  }
  if(!inter.hascov_){
    cout << name << ": plane covariance not filled" << endl;
    return -1;
  }
  // by default the covariance is not projected; the intersection must otherwise be identical
  auto nocov = interf(ktraj,surf,trange,false);
  if(nocov.hascov_ || nocov.pcov_ != PMAT() || nocov.onsurface_ != inter.onsurface_ || nocov.inbounds_ != inter.inbounds_ ||
      nocov.inrange_ != inter.inrange_ || nocov.time_ != inter.time_ || nocov.udir_ != inter.udir_ || nocov.vdir_ != inter.vdir_){
    cout << name << ": default intersection inconsistent: " << nocov << " vs " << inter << endl;
    return -1;
  }
  if(inter.inbounds_ != expectinbounds){
    cout << name << ": unexpected inbounds state " << inter.inbounds_ << endl;
    return -1;
  }
  // check the UV system is orthonormal and right-handed with the normal
  if(fabs(inter.udir_.Dot(inter.vdir_)) > ctol || fabs(inter.udir_.Dot(inter.norm_)) > ctol ||
      fabs(inter.udir_.R()-1.0) > ctol || fabs(inter.vdir_.R()-1.0) > ctol ||
      (inter.udir_.Cross(inter.vdir_) - inter.norm_).R() > ctol){
    cout << name << ": local plane UV system is not orthonormal" << endl;
    return -1;
  }
  auto const& pcov = inter.pcov_;
  if(pcov(0,0) <= 0.0 || pcov(1,1) <= 0.0 || pcov(0,0)*pcov(1,1) <= pcov(0,1)*pcov(0,1)){
    cout << name << ": plane covariance not positive definite " << pcov << endl;
    return -1;
  }
  // consistency with the trajectory planeCovariance function, using the local tangent plane
  auto tplane = surf.tangentPlane(inter.pos_);
  auto tcov = ktraj.planeCovariance(inter.time_,tplane);
  for(int irow=0;irow<2;++irow){
    for(int icol=0;icol<2;++icol){
      if(fabs(tcov(irow,icol)-pcov(irow,icol)) > ctol*sqrt(pcov(0,0)*pcov(1,1))){
        cout << name << ": inconsistent with planeCovariance " << tcov << " " << pcov << endl;
        return -1;
      }
    }
  }
  // toy MC: randomize the parameters and measure the spread of the intersection points in the local plane
  // restrict the toy intersection search to a small window around the nominal intersection, to avoid finding other
  // (physically different) intersections on other loops
  static const double dt(0.5); // ns
  TimeRange toyrange(inter.time_-dt,inter.time_+dt);
  std::mt19937 gen(seed);
  std::normal_distribution<double> gauss(0.0,1.0);
  double su(0.0), sv(0.0), suu(0.0), svv(0.0), suv(0.0);
  unsigned nsurf(0);
  for(unsigned itoy=0; itoy < ntoy; ++itoy){
    auto pars = ktraj.params();
    for(size_t ipar=0;ipar < NParams(); ++ipar)
      pars.parameters()[ipar] += gauss(gen)*sqrt(ktraj.params().covariance()(ipar,ipar));
    KTRAJ toytraj(pars,ktraj);
    auto toyinter = interf(toytraj,surf,toyrange,false); // only the position is needed
    if(!toyinter.onsurface_) continue;
    auto dpos = toyinter.pos_ - inter.pos_;
    double du = dpos.Dot(inter.udir_);
    double dv = dpos.Dot(inter.vdir_);
    su += du; sv += dv; suu += du*du; svv += dv*dv; suv += du*dv;
    ++nsurf;
  }
  if(nsurf < 0.99*ntoy){
    cout << name << ": too many toy intersections failed " << nsurf << " of " << ntoy << endl;
    return -1;
  }
  su /= nsurf; sv /= nsurf;
  double toyuu = suu/nsurf - su*su;
  double toyvv = svv/nsurf - sv*sv;
  double toyuv = suv/nsurf - su*sv;
  double toyrho = toyuv/sqrt(toyuu*toyvv);
  double rho = pcov(0,1)/sqrt(pcov(0,0)*pcov(1,1));
  if(verbose > 0){
    cout << name << " pcov  uu " << pcov(0,0) << " vv " << pcov(1,1) << " rho " << rho << endl;
    cout << name << " toyMC uu " << toyuu << " vv " << toyvv << " rho " << toyrho << " ntoy " << nsurf << endl;
  }
  if(fabs(toyuu-pcov(0,0)) > vartol*pcov(0,0) || fabs(toyvv-pcov(1,1)) > vartol*pcov(1,1) || fabs(toyrho-rho) > rhotol){
    cout << name << ": plane covariance disagrees with toy MC: pcov uu " << pcov(0,0) << " vv " << pcov(1,1) << " rho " << rho
      << " toy uu " << toyuu << " vv " << toyvv << " rho " << toyrho << endl;
    return -1;
  }
  return 0;
}

int main(int argc, char** argv) {
  unsigned ntoy(20000), seed(124223);
  int verbose(0);
  int opt;
  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"",
          long_options, &long_index )) != -1) {
    switch (opt) {
      case 'n' : ntoy = atoi(optarg);
                 break;
      case 's' : seed = atoi(optarg);
                 break;
      case 'v' : verbose = atoi(optarg);
                 break;
      default: print_usage();
               exit(EXIT_FAILURE);
    }
  }
  VEC3 bnom(0.0,0.0,1.0);
  VEC4 pos(0.0,0.0,0.0,0.0);
  MOM4 mom(40.0,30.0,60.0,0.511);
  TimeRange trange(-10.0,100.0);
  double tol(1.0e-6);
  // surfaces: tilted plane, a small rectangle on that plane (so the intersection is out of bounds), and a cylinder
  VEC3 norm = VEC3(0.3,0.2,1.0).Unit();
  VEC3 udir = norm.Cross(VEC3(0.0,1.0,0.0)).Unit();
  VEC3 pcent(5.0,-3.0,400.0);
  Plane plane(norm,udir,pcent);
  Rectangle rect(norm,udir,pcent,1.0,1.0);
  Cylinder cyl(VEC3(0.0,0.0,1.0),VEC3(0.0,0.0,0.0),200.0,1000.0);
  // trajectories with diagonal parameter covariance
  LoopHelix lhelix(pos,mom,-1,bnom,trange);
  setCovariance(lhelix,{0.5,0.5,0.5,0.5,0.002,0.05});
  CentralHelix chelix(pos,mom,-1,bnom,trange);
  setCovariance(chelix,{0.3,0.002,1.0e-5,0.3,0.001,0.05});
  KinematicLine kline(pos,mom,-1,bnom,trange);
  setCovariance(kline,{0.1,0.1,0.0005,0.0005,0.005,0.05});

  auto inter = [&](auto const& ktraj, auto const& surf, TimeRange const& range, bool pcov) { return intersect(ktraj,surf,range,tol,TimeDir::forwards,pcov); };

  int status(0);
  status |= testCovariance("LoopHelix Plane",lhelix,plane,trange,inter,true,ntoy,seed,verbose);
  status |= testCovariance("LoopHelix Rectangle (out of bounds)",lhelix,rect,trange,inter,false,ntoy,seed,verbose);
  status |= testCovariance("LoopHelix Cylinder",lhelix,cyl,trange,inter,true,ntoy,seed,verbose);
  status |= testCovariance("CentralHelix Plane",chelix,plane,trange,inter,true,ntoy,seed,verbose);
  status |= testCovariance("CentralHelix Cylinder",chelix,cyl,trange,inter,true,ntoy,seed,verbose);
  status |= testCovariance("KinematicLine Plane",kline,plane,trange,inter,true,ntoy,seed,verbose);
  status |= testCovariance("KinematicLine Rectangle (out of bounds)",kline,rect,trange,inter,false,ntoy,seed,verbose);
  status |= testCovariance("KinematicLine Cylinder",kline,cyl,trange,inter,true,ntoy,seed,verbose);
  if(status == 0)
    cout << "IntersectionCovariance test passed" << endl;
  else
    cout << "IntersectionCovariance test FAILED" << endl;
  return status;
}
