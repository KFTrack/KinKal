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
    bool first(true);
    VEC3 spos, epos;
    double tmin, tmax;
    bool startinside(true),endinside(true);
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
          TimeRange srange(std::max(trange.begin(),traj->range().begin()-gaptime),traj->range().end());
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
  Intersection<ParticleTrajectory<KinKal::KinematicLine>> intersect(ParticleTrajectory<KinKal::KinematicLine> const& kklptraj, KinKal::Surface const& surf, TimeRange trange,double tol) {
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

  template < class HELIX> Intersection<ParticleTrajectory<HELIX>> phfIntersect( ParticleTrajectory<HELIX> const& phelix, KinKal::Frustrum const& fru, TimeRange trange ,double tol) {
    // for now, call generic function.  In future, we can call the above intersection on the end disks to find the correct range more efficiently
    return pstepIntersect(phelix,fru,trange,tol);
  }

  // explicit 'specializations' for the different helix types

  Intersection<ParticleTrajectory<KinKal::LoopHelix>> intersect( ParticleTrajectory<LoopHelix> const& ploophelix, KinKal::Cylinder const& cyl, TimeRange trange ,double tol) {
    return phcIntersect(ploophelix,cyl,trange,tol);
  }
  Intersection<ParticleTrajectory<KinKal::CentralHelix>> intersect( ParticleTrajectory<CentralHelix> const& pcentralhelix, KinKal::Cylinder const& cyl, TimeRange trange ,double tol) {
    return phcIntersect(pcentralhelix,cyl,trange,tol);
  }
  Intersection<ParticleTrajectory<KinKal::LoopHelix>> intersect( ParticleTrajectory<LoopHelix> const& ploophelix, KinKal::Frustrum const& fru, TimeRange trange ,double tol) {
    return phfIntersect(ploophelix,fru,trange,tol);
  }
  Intersection<ParticleTrajectory<KinKal::CentralHelix>> intersect( ParticleTrajectory<CentralHelix> const& pcentralhelix, KinKal::Frustrum const& fru, TimeRange trange ,double tol) {
    return phfIntersect(pcentralhelix,fru,trange,tol);
  }
  Intersection<ParticleTrajectory<KinKal::LoopHelix>> intersect( ParticleTrajectory<LoopHelix> const& ploophelix, KinKal::Plane const& plane, TimeRange trange ,double tol) {
    return phpIntersect(ploophelix,plane,trange,tol);
  }
  Intersection<ParticleTrajectory<KinKal::CentralHelix>> intersect( ParticleTrajectory<CentralHelix> const& pcentralhelix, KinKal::Plane const& plane, TimeRange trange ,double tol) {
    return phpIntersect(pcentralhelix,plane,trange,tol);
  }

  // generic surface intersection; cast down till we find something that works
  template <class KTRAJ> Intersection<ParticleTrajectory<KTRAJ>> phsIntersect(ParticleTrajectory<KTRAJ> const& pktraj, Surface const& surf,TimeRange trange, double tol) {
    // use pointers to cast to avoid avoid a throw
    const Surface* surfp = &surf;
    // go through the possibilities: I don't know of anything more elegant
    auto cyl = dynamic_cast<const Cylinder*>(surfp);
    if(cyl)return intersect(pktraj,*cyl,trange,tol);
    auto fru = dynamic_cast<const Frustrum*>(surfp);
    if(fru)return intersect(pktraj,*fru,trange,tol);
    auto plane = dynamic_cast<const Plane*>(surfp);
    if(plane)return intersect(pktraj,*plane,trange,tol);
    // unknown surface subclass; return failure
    return Intersection<ParticleTrajectory<KTRAJ>>(pktraj,surf,trange,tol);
  }
  // now overload the function for helices for generic surfaces
  Intersection<ParticleTrajectory<KinKal::LoopHelix>> intersect( ParticleTrajectory<LoopHelix> const& ploophelix, KinKal::Surface const& surf, TimeRange trange ,double tol) {
    return phsIntersect(ploophelix,surf,trange,tol);
  }
  Intersection<ParticleTrajectory<KinKal::CentralHelix>> intersect( ParticleTrajectory<CentralHelix> const& pcentralhelix, KinKal::Surface const& surf, TimeRange trange ,double tol) {
    return phsIntersect(pcentralhelix,surf,trange,tol);
  }

}
#endif

