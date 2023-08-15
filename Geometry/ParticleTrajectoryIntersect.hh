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
//  Find first intersection of a particle trajectory in the specified range.  This generic implementation tests
  template <class KTRAJ, class SURF> Intersection pIntersect(ParticleTrajectory<KTRAJ> const& ptraj, SURF const& surf, TimeRange trange, double tol) {
    Intersection retval;
    // loop over pieces, and test the ones in range
    VEC3 spos, epos;
    auto curr = ptraj.nearestTraj(trange.begin());
    auto prev = curr;
    // loop until we find the best piece
    unsigned ntries(0);
    unsigned maxntries = ptraj.pieces().size(); // only try as many times as there are pieces
    do {
      ++ntries;
      // compute the intersection with the current piece
      retval = intersect(*curr,surf,trange,tol);
      if(retval.onsurface_){
        // update to use the piece nearest the current intersection time
        prev = curr;
        curr = ptraj.nearestTraj(retval.time_);
      }
    } while(curr != prev && ntries < maxntries);
    if(curr != prev){
      // we found an intersection but not on the current piece. This can happen due to gaps in the trajectory
      retval.gap_ = true;
    }
    return retval;
  }
  // KinematicLine-based particle trajectory intersect implementation can always use the generic function
  Intersection intersect(ParticleTrajectory<KinKal::KinematicLine> const& kklptraj, KinKal::Surface const& surf, TimeRange trange,double tol) {
    return pIntersect(kklptraj,surf,trange,tol);
  }

  // Helix-based particle trajectory intersect implementation with a plane
  template <class HELIX> Intersection phpIntersect(ParticleTrajectory<HELIX> const& phelix, KinKal::Plane const& plane, TimeRange trange ,double tol) {
    // for now, call generic function.  In future, we can do a smarter binary search for the correct piece using the 'constant'
    // z velocity
    return pIntersect(phelix,plane,trange,tol);
  }

  template < class HELIX> Intersection phcIntersect( ParticleTrajectory<HELIX> const& phelix, KinKal::Cylinder const& cyl, TimeRange trange ,double tol) {
    // for now, call generic function.  In future, we can call the above intersection on the end disks to find the correct range more efficiently
    return pIntersect(phelix,cyl,trange,tol);
  }

  template < class HELIX> Intersection phfIntersect( ParticleTrajectory<HELIX> const& phelix, KinKal::Frustrum const& fru, TimeRange trange ,double tol) {
    // for now, call generic function.  In future, we can call the above intersection on the end disks to find the correct range more efficiently
    return pIntersect(phelix,fru,trange,tol);
  }

  // explicit 'specializations' for the different helix types

  Intersection intersect( ParticleTrajectory<LoopHelix> const& ploophelix, KinKal::Cylinder const& cyl, TimeRange trange ,double tol) {
    return phcIntersect(ploophelix,cyl,trange,tol);
  }
  Intersection intersect( ParticleTrajectory<CentralHelix> const& pcentralhelix, KinKal::Cylinder const& cyl, TimeRange trange ,double tol) {
    return phcIntersect(pcentralhelix,cyl,trange,tol);
  }
  Intersection intersect( ParticleTrajectory<LoopHelix> const& ploophelix, KinKal::Frustrum const& fru, TimeRange trange ,double tol) {
    return phfIntersect(ploophelix,fru,trange,tol);
  }
  Intersection intersect( ParticleTrajectory<CentralHelix> const& pcentralhelix, KinKal::Frustrum const& fru, TimeRange trange ,double tol) {
    return phfIntersect(pcentralhelix,fru,trange,tol);
  }
  Intersection intersect( ParticleTrajectory<LoopHelix> const& ploophelix, KinKal::Plane const& plane, TimeRange trange ,double tol) {
    return phpIntersect(ploophelix,plane,trange,tol);
  }
  Intersection intersect( ParticleTrajectory<CentralHelix> const& pcentralhelix, KinKal::Plane const& plane, TimeRange trange ,double tol) {
    return phpIntersect(pcentralhelix,plane,trange,tol);
  }

  // generic surface intersection; cast down till we find something that works
  template <class KTRAJ> Intersection phsIntersect(ParticleTrajectory<KTRAJ> const& pktraj, Surface const& surf,TimeRange trange, double tol) {
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
    return Intersection();
  }
  // now overload the function for helices for generic surfaces
  Intersection intersect( ParticleTrajectory<LoopHelix> const& ploophelix, KinKal::Surface const& surf, TimeRange trange ,double tol) {
    return phsIntersect(ploophelix,surf,trange,tol);
  }
  Intersection intersect( ParticleTrajectory<CentralHelix> const& pcentralhelix, KinKal::Surface const& surf, TimeRange trange ,double tol) {
    return phsIntersect(pcentralhelix,surf,trange,tol);
  }

}
#endif

