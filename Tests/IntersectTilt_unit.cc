//
// Test helix intersections with cylinders and frustra when the nominal BField is tilted from the surface axis.
// Random queries are compared against a reference found by stepping over the whole time range; intersect must not
// miss a crossing the reference finds, nor return a later one.  This exercises canIntersect (KFTrack/KinKal#238).
// original author: David Brown (LBNL) 2026
//
#include "KinKal/Geometry/Intersection.hh"
#include "KinKal/Geometry/Intersect.hh"
#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Geometry/Frustrum.hh"
#include <iostream>
#include <string>
#include <random>
#include <cmath>
#include <getopt.h>

using namespace KinKal;
using std::cout;
using std::endl;

static struct option long_options[] = {
  {"nquery",     required_argument, 0, 'n'  },
  {"seed",       required_argument, 0, 's'  },
  {"verbose",    required_argument, 0, 'v'  },
  {NULL, 0,0,0}
};

void print_usage() {
  printf("Usage: IntersectTilt --nquery (number of random queries per configuration) --seed (random seed) --verbose (verbosity)\n");
}

enum SurfType {cylout=0, cylin, fruout, fruin};
static const char* surfName[4] = {"Cylinder (helix circle inside)","Cylinder (helix circle around)","Frustrum (helix circle inside)","Frustrum (helix circle around)"};

// reference: the first in-bounds crossing found by stepping, restarting past any out-of-bounds crossings
template <class HELIX, class SURF> Intersection refIntersect(HELIX const& helix, SURF const& surf, TimeRange trange, double tol, TimeDir tdir){
  static const double eps(1.0e-4); // ns: restart offset
  for(unsigned icross=0; icross < 100; ++icross){
    auto inter = stepIntersect(helix,surf,trange,tol,tdir);
    if(!inter.onsurface_ || !inter.inrange_ || inter.inbounds_) return inter;
    if(tdir == TimeDir::forwards){
      if(inter.time_+eps >= trange.end()) break;
      trange = TimeRange(inter.time_+eps,trange.end());
    } else {
      if(inter.time_-eps <= trange.begin()) break;
      trange = TimeRange(trange.begin(),inter.time_-eps);
    }
  }
  return Intersection();
}

template <class HELIX, class SURF> void compare(HELIX const& helix, SURF const& surf, TimeRange const& trange, double tol, TimeDir tdir,
    unsigned& nref, unsigned& nmiss, unsigned& nlate, unsigned& nextra, int verbose) {
  static const double dttol(1.0e-5); // ns
  auto ref = refIntersect(helix,surf,trange,tol,tdir);
  auto inter = intersect(helix,surf,trange,tol,tdir);
  if(ref.good()){
    ++nref;
    if(!inter.good()){
      ++nmiss;
      if(verbose > 1) cout << "Missed crossing: reference " << ref << endl;
    } else if((inter.time_-ref.time_)*timeDirSign(tdir) > dttol){
      ++nlate;
      if(verbose > 1) cout << "Later crossing: reference " << ref << " intersect " << inter << endl;
    }
  } else if(inter.good()) {
    ++nextra; // the reference stepping can miss grazing crossings; this is not a failure
  }
}

template <class HELIX> int testTilt(std::string const& hname, SurfType stype, double maxtilt, unsigned nquery, unsigned seed, int verbose) {
  static const double tol(1.0e-5); // mm
  std::mt19937 gen(seed);
  std::uniform_real_distribution<double> uni(0.0,1.0);
  unsigned nref(0), nmiss(0), nlate(0), nextra(0);
  for(unsigned iquery=0; iquery < nquery; ++iquery){
    // bnom tilted from z, the surface axis
    double tilt = maxtilt*uni(gen);
    double tphi = 2*M_PI*uni(gen);
    VEC3 bnom(sin(tilt)*cos(tphi),sin(tilt)*sin(tphi),cos(tilt));
    // low-momentum particle so that the helix loops within the query range
    double mom = 30.0+120.0*uni(gen);
    double cost = -0.9+1.8*uni(gen);
    double sint = sqrt(1.0-cost*cost);
    double phi = 2*M_PI*uni(gen);
    MOM4 mom4(mom*sint*cos(phi),mom*sint*sin(phi),mom*cost,0.511);
    VEC4 pos(-400.0+800.0*uni(gen),-400.0+800.0*uni(gen),-500.0+1000.0*uni(gen),0.0);
    HELIX helix(pos,mom4,-1,bnom,TimeRange(-50.0,50.0));
    // surface radius near the helix extremes, so that crossings are marginal
    double dcent = helix.center(0.0).Rho();
    double rbend = helix.bendRadius();
    double rsurf = (stype == cylin || stype == fruin) ? rbend - dcent : dcent + (uni(gen) < 0.5 ? 1.0 : -1.0)*rbend;
    rsurf += -30.0 + 60.0*uni(gen);
    double t0 = -10.0+10.0*uni(gen);
    TimeRange trange(t0,t0+2.0+18.0*uni(gen));
    auto tdir = uni(gen) < 0.5 ? TimeDir::forwards : TimeDir::backwards;
    bool isfru = stype == fruout || stype == fruin;
    double dr = -200.0 + 400.0*uni(gen); // frustrum radius change
    if(rsurf < 10.0 || (isfru && rsurf-fabs(dr) < 10.0)) continue;
    VEC3 zdir(0.0,0.0,1.0), origin(0.0,0.0,0.0);
    if(isfru){
      Frustrum fru(zdir,origin,rsurf-dr,rsurf+dr,2000.0);
      compare(helix,fru,trange,tol,tdir,nref,nmiss,nlate,nextra,verbose);
    } else {
      Cylinder cyl(zdir,origin,rsurf,2000.0);
      compare(helix,cyl,trange,tol,tdir,nref,nmiss,nlate,nextra,verbose);
    }
  }
  bool fail = nmiss > 0 || nlate > 0;
  if(verbose > 0 || fail)
    cout << hname << " " << surfName[stype] << " tilt < " << maxtilt << " rad: " << nref << " reference crossings, "
      << nmiss << " missed, " << nlate << " later, " << nextra << " extra" << (fail ? " FAILED" : "") << endl;
  return fail ? -1 : 0;
}

int main(int argc, char** argv) {
  unsigned nquery(500), seed(238);
  int verbose(0);
  int opt;
  int long_index =0;
  while ((opt = getopt_long_only(argc, argv,"",
          long_options, &long_index )) != -1) {
    switch (opt) {
      case 'n' : nquery = atoi(optarg);
                 break;
      case 's' : seed = atoi(optarg);
                 break;
      case 'v' : verbose = atoi(optarg);
                 break;
      default: print_usage();
               exit(EXIT_FAILURE);
    }
  }
  int status(0);
  for(auto stype : {cylout, cylin, fruout, fruin}){
    for(double maxtilt : {0.0, 0.05, 0.45}){
      status |= testTilt<LoopHelix>("LoopHelix",stype,maxtilt,nquery,seed,verbose);
      status |= testTilt<CentralHelix>("CentralHelix",stype,maxtilt,nquery,seed,verbose);
    }
  }
  if(status == 0)
    cout << "IntersectTilt test passed" << endl;
  else
    cout << "IntersectTilt test FAILED" << endl;
  return status;
}
