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
#include "Math/VectorUtil.h"
namespace KinKal {
  //
  // generic intersection implementation based on stepping across the surface within a given range
  //
  template <class KTRAJ, class SURF> Intersection<KTRAJ> stepIntersect( KTRAJ const& ktraj, SURF const& surf, TimeRange trange, double tol) {
    Intersection<KTRAJ> retval(ktraj,surf,trange,tol);
    double ttest = trange.begin();
    auto pos = ktraj.position3(ttest);
    double speed = ktraj.speed(ttest); // speed is constant
    bool startinside = surf.isInside(pos);
    bool stepinside;
    // set the step according to the range and tolerance.  The range division is arbitrary
    double tstep = std::max(0.05*trange.range(),tol/speed);  // trajectory range defines maximum step
    // step until we cross the surface or the time is out of range
    do {
      ttest += tstep;
      pos = ktraj.position3(ttest);
      stepinside = surf.isInside(pos);
    } while (startinside == stepinside && trange.inRange(ttest));
    if(startinside != stepinside){
      // we crossed the cylinder: backup and do a linear search.
      ttest -= tstep;
      double dist;
      do {
        retval.pos_ = ktraj.position3(ttest);
        retval.pdir_ = ktraj.direction(ttest);
        auto ray = retval.ray();
        retval.flag_ = surf.intersect(ray,dist,false,tol);
        if(retval.flag_.onsurface_){
          ttest += dist/speed;
        } else {
          break;
        }
      } while (fabs(dist) > tol);
      if(retval.flag_.onsurface_){
       // calculate the time
        retval.time_ = ttest;
        retval.pos_ = ktraj.position3(retval.time_);
        retval.pdir_ = ktraj.direction(retval.time_);
        retval.flag_.inbounds_ = surf.inBounds(retval.pos_,tol);
        retval.norm_ = surf.normal(retval.pos_);
     }
    }
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

  template < class HELIX> Intersection<HELIX> hcIntersect( HELIX const& helix, KinKal::Cylinder const& cyl, TimeRange trange ,double tol) {
    Intersection<HELIX> retval(helix,cyl,trange,tol);
    // compare directions and divide into cases
    double ddot = fabs(helix.bnom().Unit().Dot(cyl.axis()));
    if (ddot > 0.9) { // I need a more physical co-linear test TODO
      // the helix and cylinder are roughly co-linear.
      // see if the range is in bounds
      auto binb = cyl.inBounds(helix.position3(trange.begin()),tol);
      auto einb = cyl.inBounds(helix.position3(trange.end()),tol);
      if(binb && einb) {
        // both in bounds: use this range
        if(canIntersect(helix,cyl,trange,tol)){
          retval.copyResult(stepIntersect(helix,cyl,trange,tol));
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
        auto frontinter = hpIntersect(helix,frontdisk,trange,tol);
        auto backinter = hpIntersect(helix,backdisk,trange,tol);
        if(frontinter.flag_.onsurface_ && frontinter.flag_.inbounds_ && frontinter.inRange())times.push_back(frontinter.time_);
        if(backinter.flag_.onsurface_ && backinter.flag_.inbounds_ && backinter.inRange())times.push_back(backinter.time_);
        if(times.size() >=2){
          TimeRange srange(*std::min_element(times.begin(),times.end()),*std::max_element(times.begin(),times.end()));
// intersection is possible: step within the restricted range
          if(canIntersect(helix,cyl,srange,tol)){
            auto sinter = stepIntersect(helix,cyl,srange,tol);
            retval.copyResult(sinter);
          }
        }
      }
      //    } else if (ddot < 0.1) {
      //      // the helix and cylinder are mostly orthogonal, use POCA to the axis to find an initial estimate, then do a linear search
      //      // TODO. Construct a KinematicLine object from the helix axis, and a GeometricLine from the cylinder, then invoke POCA.
    } else {
    // intermediate case: use step intersection
      retval.copyResult(stepIntersect(helix,cyl,trange,tol));
    }
    return retval;
  }
  //
  // Helix and planar surfaces
  //
  template <class HELIX> Intersection<HELIX> hpIntersect( HELIX const& helix, KinKal::Plane const& plane, TimeRange trange ,double tol) {
    Intersection<HELIX> retval(helix,plane,trange,tol);
    auto axis = helix.axis(trange.begin());
    // test for the helix being circular or tangent to the plane
    double vz = helix.axisSpeed();  // speed along the helix axis
    double ddot = fabs(axis.dir_.Dot(plane.normal()));
    if(fabs(vz*trange.range()) > tol && ddot > tol ){
      // Find the intersection time of the  helix axis (along bnom) with the plane
      double dist(0.0);
      auto pinter = plane.intersect(axis,dist,true,tol);
      if(pinter.onsurface_){
        // translate the axis intersection to a time
        double tmid = trange.begin() + dist/vz;
        // bound the range of intersections by the extrema of the cylinder-plane intersection
        double tantheta = sqrt(std::max(0.0,1.0 -ddot*ddot))/ddot;
        double dt = std::max(tol/vz,helix.bendRadius()*tantheta/vz); // make range finite in case the helix is exactly co-linear with the plane normal
        // if we're already in tolerance, finish
        if(dt*vz/ddot < tol){
          retval.flag_ = pinter;
          retval.time_ = tmid;
          retval.pos_ = helix.position3(tmid);
          retval.pdir_ = helix.direction(tmid);
          retval.norm_ = plane.normal(retval.pos_);
        } else {
          TimeRange srange(tmid-dt,tmid+dt);
          if(srange.restrict(trange)){
            // step to the intersection in the restricted range.  Use a separate intersection object as the
            // range is different
            auto sinter = stepIntersect(helix,plane,srange,tol);
            retval.copyResult(sinter);
          }
        }
      }
    } else {
      // simply step to intersection
      auto sinter = stepIntersect(helix,plane,trange,tol);
      retval.copyResult(sinter);
    }
    return retval;
  }
  //
  //  Helix and frustrum
  //
  template < class HELIX> Intersection<HELIX> hfIntersect( HELIX const& helix, KinKal::Frustrum const& fru, TimeRange trange ,double tol) {
    Intersection<HELIX> retval(helix,fru,trange,tol);
    double ddot = fabs(helix.bnom().Unit().Dot(fru.axis()));
    if (ddot > 0.9) { // I need a more physical co-linear test TODO
      // the helix and frustrum are roughly co-linear.
      // see if the range is in bounds
      auto binb = fru.inBounds(helix.position3(trange.begin()),tol);
      auto einb = fru.inBounds(helix.position3(trange.end()),tol);
      if(binb && einb) {
        // both in bounds: use this range
        if(canIntersect(helix,fru,trange,tol)) {
          retval.copyResult(stepIntersect(helix,fru,trange,tol));
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
        if(frontinter.flag_.onsurface_ && frontinter.flag_.inbounds_ && frontinter.inRange())times.push_back(frontinter.time_);
        if(backinter.flag_.onsurface_ && backinter.flag_.inbounds_ && backinter.inRange())times.push_back(backinter.time_);
        if(times.size() >=2){
          TimeRange srange(*std::min_element(times.begin(),times.end()),*std::max_element(times.begin(),times.end()));
// intersection is possible: step within the restricted range
          if(canIntersect(helix,fru,srange,tol) ) {
            auto sinter = stepIntersect(helix,fru,srange,tol);
            retval.copyResult(sinter);
          }
        }
      }
      //    } else if (ddot < 0.1) {
      //      // the helix and frustrum are mostly orthogonal, use POCA to the axis to find an initial estimate, then do a linear search
      //      // TODO. Construct a KinematicLine object from the helix axis, and a GeometricLine from the frustrum, then invoke POCA.
    } else {
    // intermediate case: use step intersection
      retval.copyResult(stepIntersect(helix,fru,trange,tol));
    }
    return retval;
  }
  //
  // Tie intersect to the explicit helix implementations
  //
  Intersection<KinKal::LoopHelix> intersect( LoopHelix const& lhelix, KinKal::Cylinder const& cyl, TimeRange trange ,double tol) {
    return hcIntersect(lhelix,cyl,trange,tol);
  }
  Intersection<KinKal::CentralHelix> intersect( CentralHelix const& chelix, KinKal::Cylinder const& cyl, TimeRange trange ,double tol) {
    return hcIntersect(chelix,cyl,trange,tol);
  }
  Intersection<KinKal::LoopHelix> intersect( LoopHelix const& lhelix, KinKal::Plane const& plane, TimeRange trange ,double tol) {
    return hpIntersect(lhelix,plane,trange,tol);
  }
  Intersection<KinKal::CentralHelix> intersect( CentralHelix const& chelix, KinKal::Plane const& plane, TimeRange trange ,double tol) {
    return hpIntersect(chelix,plane,trange,tol);
  }
  Intersection<KinKal::LoopHelix> intersect( LoopHelix const& lhelix, KinKal::Frustrum const& fru, TimeRange trange ,double tol) {
    return hfIntersect(lhelix,fru,trange,tol);
  }
  Intersection<KinKal::CentralHelix> intersect( CentralHelix const& chelix, KinKal::Frustrum const& fru, TimeRange trange ,double tol) {
    return hfIntersect(chelix,fru,trange,tol);
  }
  //
  // Line trajectory with generic surfaces
  //
  template<class SURF> Intersection<KinKal::KinematicLine> intersect(KinKal::KinematicLine const& kkline, SURF const& surf, TimeRange trange,double tol) {
    Intersection<KinKal::KinematicLine> retval(kkline,surf,trange,tol);
    auto tstart = trange.begin();
    auto pos = kkline.position3(tstart);
    auto dir = kkline.direction(tstart);
    Ray ray(dir,pos);
    double dist;
    retval.flag_ = surf.intersect(ray,dist,true,tol);
    if(retval.flag_.onsurface_){
      retval.pos_ = ray.position(dist);
      retval.norm_ = surf.normal(retval.pos_);
      retval.pdir_ = dir;
      // calculate the time
      retval.time_ = tstart + dist/kkline.speed(tstart);
    }
    return retval;
  }
  // generic surface intersection; cast down till we find something that works
  // Only do this for helices, as the KinematicLine is already fully implemented for surfaces
  template <class KTRAJ> Intersection<KTRAJ> phsIntersect(KTRAJ const& ktraj, Surface const& surf,TimeRange trange, double tol) {
    // use pointers to cast to avoid avoid a throw
    const Surface* surfp = &surf;
    // go through the possibilities: I don't know of anything more elegant
    auto cyl = dynamic_cast<const Cylinder*>(surfp);
    if(cyl)return intersect(ktraj,*cyl,trange,tol);
    auto fru = dynamic_cast<const Frustrum*>(surfp);
    if(fru)return intersect(ktraj,*fru,trange,tol);
    auto plane = dynamic_cast<const Plane*>(surfp);
    if(plane)return intersect(ktraj,*plane,trange,tol);
    // unknown surface subclass; return failure
    return Intersection<KTRAJ>(ktraj,surf,trange,tol);
  }
  Intersection<KinKal::LoopHelix> intersect( LoopHelix const& ploophelix, KinKal::Surface const& surf, TimeRange trange ,double tol) {
    return phsIntersect(ploophelix,surf,trange,tol);
  }
  Intersection<KinKal::CentralHelix> intersect( CentralHelix const& pcentralhelix, KinKal::Surface const& surf, TimeRange trange ,double tol) {
    return phsIntersect(pcentralhelix,surf,trange,tol);
  }

}

#endif
