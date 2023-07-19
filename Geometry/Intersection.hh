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
#include "KinKal/Trajectory/KinematicLine.hh"

namespace KinKal {
  // intersection product
  template <class KTRAJ, class SURF> struct Intersection : public InterData {
    Intersection(KTRAJ const& ktraj, SURF const& surf,double tol) : ktraj_(ktraj), surf_(surf), tol_(tol) {}
    KTRAJ const& ktraj_; // trajectory of this intersection
    SURF const& surf_; // surf of this intersection
    double tol_; // tol used in this intersection
  };
  //
  // generic intersection implementation based on stepping across the surface within a given range
  //
  template <class KTRAJ, class SURF> Intersection<KTRAJ, SURF> stepIntersect( KTRAJ const& ktraj, SURF const& surf, TimeRange trange, double tol) {
    Intersection<KTRAJ, SURF> retval(ktraj,surf,tol);
    double ttest = trange.begin();
    auto pos = ktraj.position3(ttest);
    bool startinside = surf.isInside(pos);
    bool stepinside;
    // set the step according to curvature
    double tstep = 0.1*trange.range();  // trajectory range defines maximum step
    auto curv = surf.curvature(pos);
    if(curv > 0)tstep = std::min(tstep,0.1/(ktraj.speed()*curv));
    auto acc = ktraj.acceleration();
    if(acc > 0) tstep = std::min(tstep,0.01*ktraj.speed()/acc);
    // step until we cross the surface or the point is out-of-bounds
    do {
      ttest += tstep;
      pos = ktraj.position3(ttest);
      stepinside = surf.isInside(pos);
    } while (startinside == stepinside && surf.inBounds(pos,tol) && trange.inRange(ttest));
    if(startinside != stepinside){
      // we crossed the cylinder: backup and do a linear search.
      ttest -= tstep;
      double speed = ktraj.speed(ttest); // speed is constant
      double dist = tstep/speed;
      while (fabs(dist) > tol) {
        auto pos = ktraj.position3(ttest);
        auto dir = ktraj.direction(ttest);
        Ray ray(dir,pos);
        retval.flag_ = surf.intersect(ray,dist,false,tol);
        if(retval.flag_.onsurface_){
          ttest += dist/speed;
        } else {
          break;
        }
      }
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
  //  // specializations for different trajectory and surface types
  //  // Helix and cylinder
  //  // Implement this also for CentralHelix TODO
  //  //
  //  Intersection<KinKal::LoopHelix,KinKal::Cylinder> intersect( KinKal::LoopHelix const& lhelix, KinKal::Cylinder const& cyl, double tstart ,double tol) {
  //    Intersection<KinKal::LoopHelix, KinKal::Cylinder> retval(lhelix,cyl,tol);
  //    // compare radii and directions, and divide into cases
  //    double rratio = lhelix.rad()/cyl.radius();
  //    double ddot = fabs(lhelix.bnom().Unit().Dot(cyl.axis()));
  //    double speed = lhelix.speed(tstart); // speed is constant
  //    if (rratio < 1 && ddot*cyl.halfLength() < std::min(lhelix.rad(),cyl.radius())){
  //      // the helix is smaller than the cylinder and they are roughly co-linear.  Do a quick test to see if they could intersect within the boundary
//      // TODO
//    } else if (rratio > 1 && ddot < 0.1) {
//      // the helix and cylinder are mostly orthogonal, use POCA to the axis to find an initial estimate, then do a linear search
//      // TODO
//    } else {
//      // intermediate case: use step intersection
//      stepIntersect(lhelix,cyl,tstart,tol);
//    }
//    return retval;
//  }
//  //
//// Helix with planar surfaces can be solved generically
////
//  template <class PSURF> Intersection<KinKal::LoopHelix,class PSURF > intersect( KinKal::LoopHelix const& lhelix, PSURF const& psurf, double tstart ,double tol) {
//    Intersection<KinKal::LoopHelix, PSURF> retval(lhelix,psurf,tol);
//    // compare helix direction and plane direction, and divide into cases
//    double ddot = fabs(lhelix.bnom().Unit().Dot(cyl.axis()));
//    double speed = lhelix.speed(tstart); // speed is constant
//    if (ddot > 0.9 ){
//      // roughly colinear; use the Z component of velocity to determine an approximate time, then step.
//      // TODO
//    } else {
//      // use step intersection.  Set step according to curvature
//      double tstep = 0.01*cyl.radius()/speed;
//      stepIntersect(lhelix,psurf,tstart,tstep,tol);
//    }
//    return retval;
//  }
//


  //
  // Line trajectory can provide an exact answer for generic surfaces
  //
  template<class SURF> Intersection<KinKal::KinematicLine,SURF> intersect(KinKal::KinematicLine const& kkline, SURF const& surf, double tstart,double tol) {
    Intersection<KinKal::KinematicLine,SURF> retval(kkline,surf,tol);
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
}

#endif
