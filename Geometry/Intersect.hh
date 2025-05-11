//
//  Calculate the intersection point of a Trajectory with a surf
//  This must be specialized for every case (every pair of trajectory and surf)
//  original author: David Brown (LBNL) 2023
//
#ifndef KinKal_Intersect_hh
#define KinKal_Intersect_hh
#include "KinKal/Geometry/Intersection.hh"
#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Geometry/Frustrum.hh"
#include "KinKal/Geometry/Plane.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/CentralHelix.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/KinematicLine.hh"
#include "KinKal/General/TimeDir.hh"
#include "Math/VectorUtil.h"
#include <iostream>
namespace KinKal {
  //
  // generic intersection implementation based on stepping across the surface within a given range in the given time direction
  //
  template <class KTRAJ, class SURF> Intersection stepIntersect( KTRAJ const& ktraj, SURF const& surf, TimeRange trange, double tol,TimeDir tdir = TimeDir::forwards) {
    Intersection retval;
    double ttest = tdir == TimeDir::forwards ? trange.begin() : trange.end();
    auto pos = ktraj.position3(ttest);
    double rbend = ktraj.bendRadius();
    double tstep = trange.range();
    double speed = ktraj.speed();
    // if there's curvature, take smaller steps
    if(rbend > 0.0){
      double tlenmax = 2.0*sqrt(2.0*rbend*tol); // maximum transverse length keeping the sagitta within tolerance
      double tspeed = ktraj.transverseSpeed();
      tstep = std::min(tstep,(tlenmax+tol)/tspeed);  // trajectory range defines maximum step
    }
    // sign!
    tstep *= timeDirSign(tdir);
    bool startinside = surf.isInside(pos);
    bool stepinside;
    unsigned niter(0);
    // step until we cross the surface or the time is out of range
    do {
      ttest += tstep;
      pos = ktraj.position3(ttest);
      stepinside = surf.isInside(pos);
    } while (startinside == stepinside && trange.inRange(ttest));
    if(startinside != stepinside){
      // we crossed the surface: backup and do a linear search.
      ttest -= tstep;
      double dist = std::numeric_limits<float>::max();
      IntersectFlag rayinter;
      static const unsigned nconv(5); // when to start worrying about non-convergence
      double sumdist[2] = {0.0,0.0}; // sum of absolute distance for convergence testing
      do {
        retval.pos_ = ktraj.position3(ttest);
        retval.pdir_ = ktraj.direction(ttest);
        auto ray = retval.ray();
        rayinter = surf.intersect(ray,dist,false,tol);
        if(rayinter.onsurface_){
          ttest += dist/speed;
        } else {
          break;
        }
        ++niter;
        unsigned iconv = niter/nconv; // count convergence cycles
        if(iconv > 0){ // start testing for non-convergence
          unsigned icur = iconv % 2; // current cycle average indicator (0 or 1)
          int irem = niter - iconv*nconv;
          if( irem == 0){
            // new cycle
            if(iconv > 2){ //if we're past the 3rd cycle, test
              unsigned iprev = (iconv + 1) % 2; // previous cycle (that we just finished)
              if(sumdist[iprev]/sumdist[icur] > 0.5) {
                // claim non-convergence if the sum distance hasn't decreased by at least a factor of 2 between cycles
                // std::cout << "Non-converged step intersection ray inter dist " << dist << " avg " << sumdist[iprev] << " prev avg " << sumdist[icur] << std::endl;
                return retval; // failed
              }
            }
            sumdist[icur] = 0.0; // reset the average counter
          }
          sumdist[icur] += fabs(dist); // increatement the sum
        }
      } while (fabs(dist) > tol);
      if(rayinter.onsurface_){
        retval.onsurface_ = rayinter.onsurface_;
        retval.time_ = ttest;
        retval.pos_ = ktraj.position3(retval.time_);
        retval.pdir_ = ktraj.direction(retval.time_);
        retval.inbounds_ = surf.inBounds(retval.pos_,tol);
        retval.norm_ = surf.normal(retval.pos_);
      }
    }
    // check the final time to be in range
    retval.inrange_ = trange.inRange(retval.time_);
    return retval;
  }
  //
  // specializations for different trajectory and surface types
  // Helix and cylinder.  This must be explicit to differentiate from the line intersection below.
  // All the work is done in common with CentralHelix and LoopHelix
  //  The actual implementation is generic for any helix
  //
  //  First, a test if the circles can intersect at the specified time

  template < class HELIX> bool canIntersect( HELIX const& helix, KinKal::Cylinder const& cyl, TimeRange const& trange, double tol) {
    auto hbeg = helix.center(trange.begin());
    auto hend = helix.center(trange.end());
    // if the circles are fully contained one in the other no intersection is possible
    double dbeg = ROOT::Math::VectorUtil::Perp(hbeg - cyl.center(),cyl.axis());
    double dend = ROOT::Math::VectorUtil::Perp(hend - cyl.center(),cyl.axis());
    double drad = fabs(cyl.radius()-helix.bendRadius());
    return dbeg > drad -tol && dend > drad -tol;
  }

  template < class HELIX> bool canIntersect( HELIX const& helix, KinKal::Frustrum const& fru, TimeRange const& trange, double tol) {
    auto hbeg = helix.center(trange.begin());
    auto hend = helix.center(trange.end());
    // if the circles are fully contained one in the other no intersection is possible
    double dbeg = ROOT::Math::VectorUtil::Perp(hbeg - fru.center(),fru.axis());
    double dend = ROOT::Math::VectorUtil::Perp(hend - fru.center(),fru.axis());
    double dradbeg = fru.radius(hbeg)-helix.bendRadius();
    double dradend = fru.radius(hend)-helix.bendRadius();
    return dradbeg*dradend < 0.0 || (dbeg > fabs(dradbeg)  -tol && dend > fabs(dradend) -tol);
  }

  template < class HELIX> Intersection hcIntersect( HELIX const& helix, KinKal::Cylinder const& cyl, TimeRange trange ,double tol,TimeDir tdir = TimeDir::forwards) {
    Intersection retval;
    // test if the range is lnger than the radius. If so, make an axial approximation
    double rbend = fabs(helix.bendRadius());
    // compare directions and divide into cases
    double ddot = fabs(helix.bnom().Unit().Dot(cyl.axis()));
    if(rbend < trange.range()*helix.transverseSpeed() && ddot > 0.9) { // mostly co-linear
      // the helix and cylinder are roughly co-linear, and the helix loops within this range
      // first try to limit the range using disk-axis intersections
      // see if the range is in bounds
      auto binb = cyl.inBounds(helix.position3(trange.begin()),tol);
      auto einb = cyl.inBounds(helix.position3(trange.end()),tol);
      if(binb && einb) {
        // both in bounds: use this range
        if(canIntersect(helix,cyl,trange,tol)){
          retval = stepIntersect(helix,cyl,trange,tol,tdir);
        }
      } else {
        // create a list of times to select the most restrictive range
        std::vector<double>times;
        if(binb)times.push_back(trange.begin());
        if(einb)times.push_back(trange.end());
        // find the intersections with the front and back disks.  Use those to refine the range
        // note we don't know the orientation of the trajectory WRT the axis so we have to try both
        auto frontdisk = cyl.frontDisk();
        auto backdisk = cyl.backDisk();
        auto frontinter = hpIntersect(helix,frontdisk,trange,tol,tdir);
        auto backinter = hpIntersect(helix,backdisk,trange,tol,tdir);
        if(trange.inRange(frontinter.time_))times.push_back(frontinter.time_);
        if(trange.inRange(backinter.time_))times.push_back(backinter.time_);
        if(times.size() >=2){
          TimeRange srange(*std::min_element(times.begin(),times.end()),*std::max_element(times.begin(),times.end()));
          // intersection is possible: step within the restricted range
          if(canIntersect(helix,cyl,srange,tol)){
            retval = stepIntersect(helix,cyl,srange,tol,tdir);
          }
        }
      }
    } else {
      // intermediate case: use step intersection
      retval = stepIntersect(helix,cyl,trange,tol,tdir);
    }
    return retval;
  }
  //
  // Helix and planar surfaces
  //
  template <class HELIX> Intersection hpIntersect( HELIX const& helix, KinKal::Plane const& plane, TimeRange trange ,double tol,TimeDir tdir = TimeDir::forwards) {
    Intersection retval;
    // test if the range is lnger than the radius. If so, make an axial approximation
    double rbend = fabs(helix.bendRadius());
    if(rbend > trange.range()*helix.transverseSpeed()){
      retval = stepIntersect(helix,plane,trange,tol,tdir);
    } else {
      // the trajectory curves macroscopically over this range: see if the axis is normal
      double tstart = tdir == TimeDir::forwards ? trange.begin() : trange.end();
      auto ray = helix.linearize(tstart);
      auto velo = helix.velocity(tstart);
      if(tdir == TimeDir::backwards) ray.reverse(); // reverse if going backwards in time
      double vax = velo.Dot(ray.direction()); // speed along axis
      // test for the helix being circular or tangent to the plane
      double ddot = fabs(ray.direction().Dot(plane.normal()));
      double zrange = fabs(vax*trange.range());
      if(zrange > tol && ddot > tol/zrange ){
        // Find the intersection time of the  helix ray (along bnom) with the plane to reduce the range
        double dist(0.0);
        auto pinter = plane.intersect(ray,dist,true,tol);
        if(pinter.onsurface_){
          // translate the ray intersection distance to a time
          double tmid = tstart + dist/vax;
          // bound the range of intersections by the extrema of the cylinder-plane intersection
          double tantheta = sqrt(std::max(1.0 -ddot*ddot,0.0))/ddot;
          // distance along axis to the surface to bound the reduced range; add tolerance
          double dd = tol+rbend*tantheta;
          // intersection should be bounded by the
          auto dt = fabs(dd/vax);
          TimeRange srange(tmid-dt,tmid+dt);
          // if this range overlaps with the original, compute the intersection
          if(srange.restrict(trange)){
            // step to the intersection in the restricted range.  Use a separate intersection object as the
            // range is different
            retval = stepIntersect(helix,plane,srange,tol,tdir);
          }
        }
      } else {
        // simply step to intersection
        retval = stepIntersect(helix,plane,trange,tol,tdir);
      }
    }
    return retval;
  }
  //
  //  Helix and frustrum
  //
  template < class HELIX> Intersection hfIntersect( HELIX const& helix, KinKal::Frustrum const& fru, TimeRange trange ,double tol,TimeDir tdir = TimeDir::forwards) {
    Intersection retval;
    double ddot = fabs(helix.bnom().Unit().Dot(fru.axis()));
    if (ddot > 0.9) { // I need a more physical co-linear test TODO
      // the helix and frustrum are roughly co-linear.
      // see if the range is in bounds
      auto binb = fru.inBounds(helix.position3(trange.begin()),tol);
      auto einb = fru.inBounds(helix.position3(trange.end()),tol);
      if(binb && einb) {
        // both in bounds: use this range
        if(canIntersect(helix,fru,trange,tol)) {
          retval  = stepIntersect(helix,fru,trange,tol,tdir);
        }
      } else {
        // create a list of times to select the most restrictive range
        std::vector<double>times;
        if(binb)times.push_back(trange.begin());
        if(einb)times.push_back(trange.end());
        // find the intersections with the front and back disks.  Use those to refine the range
        // note we don't know the orientation of the trajectory WRT the axis so we have to try both
        auto frontdisk = fru.frontDisk();
        auto backdisk = fru.backDisk();
        auto frontinter = hpIntersect(helix,frontdisk,trange,tol);
        auto backinter = hpIntersect(helix,backdisk,trange,tol);
        if(trange.inRange(frontinter.time_))times.push_back(frontinter.time_);
        if(trange.inRange(backinter.time_))times.push_back(backinter.time_);
        if(times.size() >=2){
          TimeRange srange(*std::min_element(times.begin(),times.end()),*std::max_element(times.begin(),times.end()));
          // intersection is possible: step within the restricted range
          if(canIntersect(helix,fru,srange,tol) ) {
            retval = stepIntersect(helix,fru,srange,tol,tdir);
          }
        }
      }
      //    } else if (ddot < 0.1) {
      //      // the helix and frustrum are mostly orthogonal, use POCA to the axis to find an initial estimate, then do a linear search
      //      // TODO. Construct a KinematicLine object from the helix axis, and a GeometricLine from the frustrum, then invoke POCA.
  } else {
    // intermediate case: use step intersection
    retval  = stepIntersect(helix,fru,trange,tol,tdir);
  }
  return retval;
  }
  //
  // Tie intersect to the explicit helix implementations
  //
  Intersection intersect( LoopHelix const& lhelix, KinKal::Cylinder const& cyl, TimeRange trange ,double tol,TimeDir tdir = TimeDir::forwards) {
    return hcIntersect(lhelix,cyl,trange,tol,tdir);
  }
  Intersection intersect( CentralHelix const& chelix, KinKal::Cylinder const& cyl, TimeRange trange ,double tol,TimeDir tdir = TimeDir::forwards) {
    return hcIntersect(chelix,cyl,trange,tol,tdir);
  }
  Intersection intersect( LoopHelix const& lhelix, KinKal::Plane const& plane, TimeRange trange ,double tol,TimeDir tdir = TimeDir::forwards) {
    return hpIntersect(lhelix,plane,trange,tol,tdir);
  }
  Intersection intersect( CentralHelix const& chelix, KinKal::Plane const& plane, TimeRange trange ,double tol,TimeDir tdir = TimeDir::forwards) {
    return hpIntersect(chelix,plane,trange,tol,tdir);
  }
  Intersection intersect( LoopHelix const& lhelix, KinKal::Frustrum const& fru, TimeRange trange ,double tol,TimeDir tdir = TimeDir::forwards) {
    return hfIntersect(lhelix,fru,trange,tol,tdir);
  }
  Intersection intersect( CentralHelix const& chelix, KinKal::Frustrum const& fru, TimeRange trange ,double tol,TimeDir tdir = TimeDir::forwards) {
    return hfIntersect(chelix,fru,trange,tol,tdir);
  }
  //
  // Line trajectory with generic surfaces; no specialization needed since all surfaces provide a ray intersection implementation
  //
  Intersection intersect(KinKal::KinematicLine const& kkline, KinKal::Surface const& surf, TimeRange trange,double tol,TimeDir tdir = TimeDir::forwards) {
    Intersection retval;
    auto tstart = tdir == TimeDir::forwards ? trange.begin() : trange.end();
    auto pos = kkline.position3(tstart);
    auto dir = kkline.direction(tstart)*timeDirSign(tdir);
    Ray ray(dir,pos);
    double dist;
    auto rayinter = surf.intersect(ray,dist,true,tol);
    if(rayinter.onsurface_){
      retval.onsurface_ = rayinter.onsurface_;
      retval.inbounds_ = rayinter.inbounds_;
      retval.pos_ = ray.position(dist);
      retval.norm_ = surf.normal(retval.pos_);
      retval.pdir_ = dir;
      // calculate the time
      retval.time_ = tstart + dist*timeDirSign(tdir)/kkline.speed(tstart);
    }
    // check the final time to be in range;
    retval.inrange_ = trange.inRange(retval.time_);
    return retval;
  }
  // generic surface intersection cast down till we find something that works.  This will only be used for helices, as KinematicLine
  // is already complete
  template <class KTRAJ> Intersection hsIntersect(KTRAJ const& ktraj, Surface const& surf,TimeRange trange, double tol,TimeDir tdir = TimeDir::forwards) {
    // use pointers to cast to avoid avoid a throw
    const Surface* surfp = &surf;
    // go through the possibilities: I don't know of anything more elegant
    auto plane = dynamic_cast<const Plane*>(surfp);
    if(plane)return intersect(ktraj,*plane,trange,tol,tdir);
    auto cyl = dynamic_cast<const Cylinder*>(surfp);
    if(cyl)return intersect(ktraj,*cyl,trange,tol,tdir);
    auto fru = dynamic_cast<const Frustrum*>(surfp);
    if(fru)return intersect(ktraj,*fru,trange,tol,tdir);
    // unknown surface subclass; return failure
    return Intersection();
  }
  // now provide the explicit generic interface
  Intersection intersect( LoopHelix const& ploophelix, KinKal::Surface const& surf, TimeRange trange ,double tol,TimeDir tdir = TimeDir::forwards) {
    return hsIntersect(ploophelix,surf,trange,tol,tdir);
  }
  Intersection intersect( CentralHelix const& pcentralhelix, KinKal::Surface const& surf, TimeRange trange ,double tol,TimeDir tdir = TimeDir::forwards) {
    return hsIntersect(pcentralhelix,surf,trange,tol,tdir);
  }

}

#endif
