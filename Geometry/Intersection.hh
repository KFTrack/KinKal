//
//  Calculate the intersection point of a Trajectory with a surf
//  This must be specialized for every case (every pair of trajectory and surf)
//  original author: David Brown (LBNL) 2023
//
#ifndef KinKal_Intersection_hh
#define KinKal_Intersection_hh
#include "KinKal/Geometry/InterData.hh"
#include "KinKal/Geometry/Cylinder.hh"
#include "KinKal/Geometry/Plane.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/CentralHelix.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/KinematicLine.hh"
#include "Math/VectorUtil.h"

namespace KinKal {
  // intersection product
  template <class KTRAJ> struct Intersection : public InterData {
    Intersection(KTRAJ const& ktraj, Surface const& surf,TimeRange const& trange, double tol) : InterData(trange), ktraj_(ktraj), surf_(surf), tol_(tol) {}
    KTRAJ const& ktraj_; // trajectory of this intersection
    Surface const& surf_; // surf of this intersection
    double tol_; // tol used in this intersection
  };
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
  template < class HELIX> Intersection<HELIX> helixIntersect( HELIX const& helix, KinKal::Cylinder const& cyl, TimeRange trange ,double tol) {
    // compare directions and divide into cases
    double ddot = fabs(helix.bnom().Unit().Dot(cyl.axis()));
    if (ddot > 0.9) {
      // the helix and cylinder are roughly co-linear.
      // find the intersections with the front and back disks.  Use that to refine the range
      auto frontdisk = cyl.frontDisk();
      auto backdisk = cyl.backDisk();
      auto frontinter = planeIntersect(helix,frontdisk,trange,tol);
      auto backinter = planeIntersect(helix,frontdisk,trange,tol);
      if(frontinter.flag_.onsurface_ && backinter.flag_.onsurface_){
        // front and back disks intersected.  Use these to define a restricted range.
        double tmin = std::min(frontinter.time_,backinter.time_);
        double tmax = std::max(frontinter.time_,backinter.time_);
        TimeRange srange(std::max(tmin,trange.begin()),std::min(tmax,trange.end()));
        // do a rough test if an intersection is possible by comparing the positions at the extrema
        auto axis = helix.axis(srange.begin());
        auto vz = helix.velocity(srange.begin()).Dot(axis.dir_);
        double dist = srange.range()*vz;
        auto backpos = axis.position(dist);
        // if the circles at both ends are fully contained one in the other no intersection is possible
        double drfront = ROOT::Math::VectorUtil::Perp(axis.start_ - cyl.center(),cyl.axis());
        double drback = ROOT::Math::VectorUtil::Perp(backpos - cyl.center(),cyl.axis());
        double dr = fabs(cyl.radius()-helix.bendRadius());
        if(drfront < dr && drback < dr){
          // intersection is possible: step within the restricted range
          return stepIntersect(helix,cyl,srange,tol);
        }
      }
//    } else if (ddot < 0.1) {
//      // the helix and cylinder are mostly orthogonal, use POCA to the axis to find an initial estimate, then do a linear search
//      // TODO. Construct a KinematicLine object from the helix axis, and a GeometricLine from the cylinder, then invoke POCA.
    } else {
      // intermediate case: use step intersection
      return stepIntersect(helix,cyl,trange,tol);
    }
    return Intersection<HELIX> (helix,cyl,trange,tol);
  }
  //
  // Helix and planar surfaces
  //
  template <class HELIX> Intersection<HELIX> planeIntersect( HELIX const& helix, KinKal::Plane const& plane, TimeRange trange ,double tol) {
    Intersection<HELIX> retval(helix,plane,trange,tol);
   // Find the intersection time of the  helix axis (along bnom) with the plane
    auto axis = helix.axis(trange.begin());
    double ddot = fabs(axis.dir_.Dot(plane.normal()));
    double vz = helix.velocity(trange.mid()).Dot(axis.dir_);
    double dist(0.0);
    auto pinter = plane.intersect(axis,dist,true,tol);
    if(pinter.onsurface_){
      double tmid = trange.begin() + dist/vz;
    // use the difference in axis direction and curvature to bound the range of intersection times
      double tantheta = sqrt(std::max(0.0,1.0 -ddot*ddot))/ddot;
      double dt = std::max(tol,helix.bendRadius()*tantheta)/vz; // make range finite in case the helix is exactly co-linear with the plane normal
      TimeRange srange(tmid-dt,tmid+dt);
      // now step to the exact intersection
      auto sinter = stepIntersect(helix,plane,srange,tol);
      retval.time_ = sinter.time_;
      retval.pos_  = sinter.pos_;
      retval.pdir_ = sinter.pdir_;
      retval.flag_ = sinter.flag_;
      retval.norm_ = sinter.norm_;
    }
    return retval;
  }
  //
  Intersection<KinKal::LoopHelix> intersect( LoopHelix const& lhelix, KinKal::Cylinder const& cyl, TimeRange trange ,double tol) {
    return helixIntersect(lhelix,cyl,trange,tol);
  }
  Intersection<KinKal::CentralHelix> intersect( CentralHelix const& chelix, KinKal::Cylinder const& cyl, TimeRange trange ,double tol) {
    return helixIntersect(chelix,cyl,trange,tol);
  }
  Intersection<KinKal::LoopHelix> intersect( LoopHelix const& lplane, KinKal::Plane const& plane, TimeRange trange ,double tol) {
    return planeIntersect(lplane,plane,trange,tol);
  }
  Intersection<KinKal::CentralHelix> intersect( CentralHelix const& cplane, KinKal::Plane const& plane, TimeRange trange ,double tol) {
    return planeIntersect(cplane,plane,trange,tol);
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

//  Find first intersection of a particle trajectory.  First, generic implementation looping over pieces
  template <class KTRAJ, class SURF> Intersection<ParticleTrajectory<KTRAJ>> pieceIntersect(ParticleTrajectory<KTRAJ> const& ptraj, SURF const& surf, TimeRange trange, double tol) {
    Intersection<ParticleTrajectory<KTRAJ>> retval(ptraj,surf,trange,tol);
    // loop over pieces, and test the ones in range
    bool first(false);
    bool startinside,endinside;
    for(auto traj : ptraj.pieces()) {
      if(trange.inRange(traj->range().begin()) || trange.inRange(traj->range().end())){
        if(first){
          double tmin = std::max(trange.begin(),traj->range().begin());
          auto spos = traj->position3(tmin);
          startinside = surf.isInside(spos);
          first = false;
        }
        double tmax = std::min(trange.end(),traj->range().end());
        auto epos = traj->position3(tmax);
        endinside = surf.isInside(epos);
        if(startinside != endinside){
          // we crossed the surface: find the exact intersection
          TimeRange srange(std::max(trange.begin(),traj->range().begin()),tmax);
          auto pinter = intersect(*traj,surf,srange,tol);
          if(pinter.flag_.onsurface_ && pinter.inRange()){
            retval.flag_ = pinter.flag_;
            retval.pos_ = pinter.pos_;
            retval.norm_ = pinter.norm_;
            retval.pdir_ = pinter.pdir_;
            retval.time_ = pinter.time_;
          }
        }
      }
    }
    return retval;
  }

}
//    auto curv = surf.curvature(pos);
//    if(curv > 0)tstep = std::min(tstep,0.1/(ktraj.speed()*curv));
//    auto acc = ktraj.acceleration();
//    if(acc > 0) tstep = std::min(tstep,0.01*ktraj.speed()/acc);

#endif
