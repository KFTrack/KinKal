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
  //  Find first intersection of a particle trajectory in the specified range.  This is a generic implementation, dependent on a
  //  reasonable estimate of the time.
  template <class KTRAJ, class SURF> Intersection pIntersect(ParticleTrajectory<KTRAJ> const& ptraj, SURF const& surf, TimeRange trange, double tstart, double tol,TimeDir tdir = TimeDir::forwards) {
    Intersection retval;
    // loop over pieces, and test the ones in range
    VEC3 spos, epos;
    auto curr = ptraj.nearestTraj(tstart);
    auto prev = curr;
    // loop until we find the best piece
    unsigned ntries(0);
    unsigned maxntries = ptraj.pieces().size(); // only try as many times as there are pieces
    do {
      ++ntries;
      // compute the intersection with the current piece
      retval = intersect(*curr,surf,trange,tol,tdir);
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

  // estimate starting time for ptraj-plane intersections. This helps localize which piece to use.
  template <class KTRAJ> double startTime(ParticleTrajectory<KTRAJ> const& pktraj, KinKal::Plane const& plane, TimeRange trange ,double tol,TimeDir tdir) {
    // initialize starting time. This is a bootstrap, no precision needed
    double tstart = (tdir == TimeDir::forwards) ? trange.begin() : trange.end();
    auto ttraj = pktraj.nearestTraj(tstart);
    // linear approximation to the trajectory at that end
    auto ray = ttraj->linearize(tstart);
    if(tdir == TimeDir::backwards)ray.reverse();
    double dist; // signed distance from ray start to the plane
    auto ainter = plane.intersect(ray,dist,false,tol);
    if(ainter.onsurface_){
      auto velo = ttraj->velocity(tstart);
      double vax = velo.Dot(ray.direction());
      tstart += dist/vax;
    }
    return tstart;
  }

  template <class KTRAJ, class SURF> KinKal::Disk startDisk(ParticleTrajectory<KTRAJ> const& ptraj, SURF const& surf, TimeRange trange, double tol,TimeDir tdir) {
    auto fplane = surf.frontDisk();
    auto bplane = surf.backDisk();
    double tstart = (tdir == TimeDir::forwards) ? trange.begin() : trange.end();
    KinKal::Ray ray = (tdir == TimeDir::forwards) ?  ptraj.front().linearize(tstart) : ptraj.back().linearize(tstart);
    if(tdir != TimeDir::forwards)ray.reverse();
    auto fdist = (ray.start() - fplane.center()).Dot(ray.direction());
    auto bdist = (ray.start() - bplane.center()).Dot(ray.direction());
    // choose the closest positive
    if(fdist < bdist && fdist > 0.0)
      return fplane;
    else if(bdist < fdist && bdist > 0.0)
      return bplane;
    else
      return surf.midDisk();
  }

  // KinematicLine-based particle trajectory intersect implementation can always use the generic function
  Intersection intersect(ParticleTrajectory<KinKal::KinematicLine> const& kklptraj, KinKal::Surface const& surf, TimeRange trange, double tol,TimeDir tdir = TimeDir::forwards) {
    return pIntersect(kklptraj,surf,trange,trange.begin(),tol,tdir);
  }

  // Helix-based particle trajectory intersect implementation with a plane
  template <class HELIX> Intersection phpIntersect(ParticleTrajectory<HELIX> const& phelix, KinKal::Plane const& plane, TimeRange trange ,double tol,TimeDir tdir = TimeDir::forwards) {
    auto tstart = startTime(phelix,plane,trange,tol,tdir);
    return pIntersect(phelix,plane,trange,tstart,tol,tdir);
  }

  template < class HELIX> Intersection phcIntersect( ParticleTrajectory<HELIX> const& phelix, KinKal::Cylinder const& cyl, TimeRange trange ,double tol, TimeDir tdir = TimeDir::forwards) {
    auto disk = startDisk(phelix,cyl,trange,tol,tdir);
    double tstart = startTime(phelix,disk,trange,tol,tdir);
    return pIntersect(phelix,cyl,trange,tstart,tol,tdir);
  }

  template < class HELIX> Intersection phfIntersect( ParticleTrajectory<HELIX> const& phelix, KinKal::Frustrum const& fru, TimeRange trange ,double tol, TimeDir tdir = TimeDir::forwards) {
    auto disk = startDisk(phelix,fru,trange,tol,tdir);
    double tstart = startTime(phelix,disk,trange,tol,tdir);
    return pIntersect(phelix,fru,trange,tstart,tol,tdir);
  }
  // explicit 'specializations' for the different helix types

  Intersection intersect( ParticleTrajectory<LoopHelix> const& ploophelix, KinKal::Cylinder const& cyl, TimeRange trange ,double tol, TimeDir tdir = TimeDir::forwards) {
    return phcIntersect(ploophelix,cyl,trange,tol,tdir);
  }
  Intersection intersect( ParticleTrajectory<CentralHelix> const& pcentralhelix, KinKal::Cylinder const& cyl, TimeRange trange ,double tol, TimeDir tdir = TimeDir::forwards) {
    return phcIntersect(pcentralhelix,cyl,trange,tol,tdir);
  }
  Intersection intersect( ParticleTrajectory<LoopHelix> const& ploophelix, KinKal::Frustrum const& fru, TimeRange trange ,double tol, TimeDir tdir = TimeDir::forwards) {
    return phfIntersect(ploophelix,fru,trange,tol,tdir);
  }
  Intersection intersect( ParticleTrajectory<CentralHelix> const& pcentralhelix, KinKal::Frustrum const& fru, TimeRange trange ,double tol, TimeDir tdir = TimeDir::forwards) {
    return phfIntersect(pcentralhelix,fru,trange,tol,tdir);
  }
  Intersection intersect( ParticleTrajectory<LoopHelix> const& ploophelix, KinKal::Plane const& plane, TimeRange trange ,double tol, TimeDir tdir = TimeDir::forwards) {
    return phpIntersect(ploophelix,plane,trange,tol,tdir);
  }
  Intersection intersect( ParticleTrajectory<CentralHelix> const& pcentralhelix, KinKal::Plane const& plane, TimeRange trange ,double tol, TimeDir tdir = TimeDir::forwards) {
    return phpIntersect(pcentralhelix,plane,trange,tol,tdir);
  }

  // generic surface intersection; cast down till we find something that works
  template <class KTRAJ> Intersection phsIntersect(ParticleTrajectory<KTRAJ> const& pktraj, Surface const& surf,TimeRange trange, double tol, TimeDir tdir = TimeDir::forwards) {
    // use pointers to cast to avoid avoid a throw
    const Surface* surfp = &surf;
    // go through the possibilities: I don't know of anything more elegant
    auto plane = dynamic_cast<const Plane*>(surfp);
    if(plane)return intersect(pktraj,*plane,trange,tol,tdir);
    auto cyl = dynamic_cast<const Cylinder*>(surfp);
    if(cyl)return intersect(pktraj,*cyl,trange,tol,tdir);
    auto fru = dynamic_cast<const Frustrum*>(surfp);
    if(fru)return intersect(pktraj,*fru,trange,tol,tdir);
    // unknown surface subclass; return failure
    return Intersection();
  }
  // now overload the function for helices for generic surfaces
  Intersection intersect( ParticleTrajectory<LoopHelix> const& ploophelix, KinKal::Surface const& surf, TimeRange trange ,double tol, TimeDir tdir = TimeDir::forwards) {
    return phsIntersect(ploophelix,surf,trange,tol,tdir);
  }
  Intersection intersect( ParticleTrajectory<CentralHelix> const& pcentralhelix, KinKal::Surface const& surf, TimeRange trange ,double tol, TimeDir tdir = TimeDir::forwards) {
    return phsIntersect(pcentralhelix,surf,trange,tol,tdir);
  }

}
#endif

