//
//  Calculate the intersection point of a ParticleTrajectory with a surface
//  This must be specialized for every case (every pair of trajectory and surf)
//  original author: David Brown (LBNL) 2023
//
#ifndef KinKal_ParticleTrajectoryIntersect_hh
#define KinKal_ParticleTrajectoryIntersect_hh
#include "KinKal/Geometry/Intersection.hh"
#include "KinKal/Geometry/Intersect.hh"
namespace KinKal {
  // find intersection by stepping over pieces. This works when the pieces have small curvature
  template <class KTRAJ, class SURF> Intersection intersect(ParticleTrajectory<KTRAJ> const& ptraj, SURF const& surf, TimeRange trange, double tol,TimeDir tdir = TimeDir::forwards) {
    Intersection retval;
    double tstart = (tdir == TimeDir::forwards) ? trange.begin() : trange.end();
    double tend = (tdir == TimeDir::forwards) ? trange.end() : trange.begin();
    size_t istart = ptraj.nearestIndex(tstart);
    size_t iend = ptraj.nearestIndex(tend);
    if(istart == iend){
      retval = intersect(ptraj.piece(istart),surf,trange,tol,tdir);
    } else {
      int istep = (tdir == TimeDir::forwards) ? 1 : -1;
      do {
        // if we can approximate this piece as a line, simply test at the endpoints. Time order doesn't matter
        bool testinter(true);
        auto ttraj = ptraj.indexTraj(istart);
        double sag = ttraj->sagitta(ttraj->range().range());
        if(sag < tol){
          auto spos = ttraj->position3(ttraj->range().begin());
          bool sinside = surf.isInside(spos);
          auto epos = ttraj->position3(ttraj->range().end());
          bool einside = surf.isInside(epos);
          testinter = sinside != einside;
        }
        if(testinter){
          // try to intersect. use a temporary
          auto srange = ttraj->range();
          if(srange.restrict(trange)){
            auto tinter = intersect(*ttraj,surf,srange,tol,tdir);
            if(tinter.good()){
              retval = tinter;
              break;
            }
          }
        }
        // update
        istart += istep;
      } while( istart != iend+istep); // include the last piece in the test
    }
    return retval;
  }
}
#endif

