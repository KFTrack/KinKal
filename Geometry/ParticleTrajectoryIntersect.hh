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
//  Find first intersection of a particle trajectory.  This is a generic implementation looping over pieces
//  It can miss double-intersections on the same piece; those are special cases that need to be tested for
//  in a dedicated function
  template <class KTRAJ, class SURF> Intersection<ParticleTrajectory<KTRAJ>> pstepIntersect(ParticleTrajectory<KTRAJ> const& ptraj, SURF const& surf, TimeRange trange, double tol) {
    Intersection<ParticleTrajectory<KTRAJ>> retval(ptraj,surf,trange,tol);
    // loop over pieces, and test the ones in range
    bool first(false);
    VEC3 spos, epos;
    double tmin, tmax;
    bool startinside,endinside;
    for(auto traj : ptraj.pieces()) {
      if(trange.inRange(traj->range().begin()) || trange.inRange(traj->range().end())){
        tmin = std::max(trange.begin(),traj->range().begin());
        spos = traj->position3(tmin);
        startinside = surf.isInside(spos);
        if(first){ // preload first valid end side
          first = false;
          endinside = startinside;
        }
        // compare to the previous end position
        if(startinside != endinside){
          // we crossed the surface through a trajectory gap. Search for an intersection on this piece
          // including an early buffer for the gap
          double gaptime = (spos-epos).R()/traj->speed();
          TimeRange srange(std::max(trange.begin(),traj->range().begin()-gaptime),tmax);
          auto pinter = intersect(*traj,surf,srange,tol);
          if(pinter.flag_.onsurface_ && pinter.inRange()){
            // we found the intersection; set return value and finish
            retval.copyResult(pinter);
            break;
          }
        }
        tmax = std::min(trange.end(),traj->range().end());
        epos = traj->position3(tmax);
        endinside = surf.isInside(epos);
        // test for crossing in this piece
        if(startinside != endinside){
          // we crossed the surface: find the exact intersection
          TimeRange srange(std::max(trange.begin(),traj->range().begin()),tmax);
          auto pinter = intersect(*traj,surf,srange,tol);
          if(pinter.flag_.onsurface_ && pinter.inRange()){
            // we found the intersection; set return value and finish
            retval.copyResult(pinter);
            break;
          }
        }
      }
    }
    return retval;
  }
  // KinematicLine-based particle trajectory intersect implementation can always use the generic function
  template<class SURF> Intersection<ParticleTrajectory<KinKal::KinematicLine>> intersect(ParticleTrajectory<KinKal::KinematicLine> const& kklptraj, SURF const& surf, TimeRange trange,double tol) {
    return pstepIntersect(kklptraj,surf,trange,tol);
  }

  // Helix-based particle trajectory intersect implementation with a plane
  template <class HELIX> Intersection<ParticleTrajectory<HELIX>> phpIntersect(ParticleTrajectory<HELIX> const& phelix, KinKal::Plane const& plane, TimeRange trange ,double tol) {
    // for now, call generic function.  In future, we can do a smarter binary search for the correct piece using the 'constant'
    // z velocity
    return pstepIntersect(phelix,plane,trange,tol);
  }

  template < class HELIX> Intersection<ParticleTrajectory<HELIX>> phcIntersect( ParticleTrajectory<HELIX> const& phelix, KinKal::Cylinder const& cyl, TimeRange trange ,double tol) {
    // for now, call generic function.  In future, we can call the above intersection on the end disks to find the correct range more efficiently
    return pstepIntersect(phelix,cyl,trange,tol);
  }

  // explicit 'specializations' for the different helix types

  Intersection<ParticleTrajectory<KinKal::LoopHelix>> intersect( ParticleTrajectory<LoopHelix> const& ploophelix, KinKal::Cylinder const& cyl, TimeRange trange ,double tol) {
    return phcIntersect(ploophelix,cyl,trange,tol);
  }
  Intersection<ParticleTrajectory<KinKal::CentralHelix>> intersect( ParticleTrajectory<CentralHelix> const& pcentralhelix, KinKal::Cylinder const& cyl, TimeRange trange ,double tol) {
    return phcIntersect(pcentralhelix,cyl,trange,tol);
  }
  Intersection<ParticleTrajectory<KinKal::LoopHelix>> intersect( ParticleTrajectory<LoopHelix> const& ploophelix, KinKal::Plane const& plane, TimeRange trange ,double tol) {
    return phpIntersect(ploophelix,plane,trange,tol);
  }
  Intersection<ParticleTrajectory<KinKal::CentralHelix>> intersect( ParticleTrajectory<CentralHelix> const& pcentralhelix, KinKal::Plane const& plane, TimeRange trange ,double tol) {
    return phpIntersect(pcentralhelix,plane,trange,tol);
  }


}
#endif
