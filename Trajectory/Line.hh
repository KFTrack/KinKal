#ifndef KinKal_Line_hh
#define KinKal_Line_hh

//
//  Linear time-based trajectory with a constant velocity.
//  Used as part of the kinematic Kalman fit
//
#include <memory>
#include "KinKal/General/Vectors.hh"
#include "KinKal/General/TimeRange.hh"
#include "KinKal/Trajectory/DistanceToTime.hh"
#include "KinKal/Trajectory/GeometricLine.hh"

namespace KinKal {
  class Line {
    public:
      // construct from a spacetime point (typically the measurement position and time) and propagation velocity (mm/ns).
      Line(VEC4 const& p0, VEC3 const& svel, double length);
      Line(VEC3 const& p0, double t0, VEC3 const& svel, double length);
      // construct from 2 points plus timing information.  P0 is the measurement (near) end, p1 the far end.  Signals propagate from far to near
      Line(VEC3 const& p0, VEC3 const& p1, double t0, double speed );
      Line(VEC3 const& p0, double length, VEC3 const& svel, double t0, std::shared_ptr<DistanceToTime> d2t);
      // accessors
      double t0() const { return t0_; }
      double& t0() { return t0_; } // detector updates need to refine t0
      // signal ends at pos0
      VEC3 const& startPosition() const { return gline_.startPosition(); }
      VEC3 endPosition() const { return gline_.endPosition() ; }
      double speed(double time) const { return d2t_->speed(d2t_->distance(time-t0_)); }
      double length() const { return gline_.length(); }
      VEC3 const& direction() const { return gline_.direction(); }
      // TOCA to a point
      double TOCA(VEC3 const& point) const;
      // geometric accessors
      VEC3 position3(double time) const;
      VEC4 position4(double time) const;
      VEC3 velocity(double time) const;
      VEC3 const& direction(double time) const { return gline_.direction(); }
      void print(std::ostream& ost, int detail) const;
      double timeAtMidpoint() const { return t0_ + d2t_->time(0.5*length()); }

    private:
      //VEC3 pos0_, dir_; // position and direction
       double t0_; // intial time (at pos0)
       //double speed_; // signed linear velocity, translates time to distance along the trajectory (mm/nsec)
      //double length_; // line length
      std::shared_ptr<DistanceToTime> d2t_; // represents the possibly nonlinear distance to time relationship of the line
      GeometricLine gline_; // geometic representation of the line
  };
  std::ostream& operator <<(std::ostream& ost, Line const& tline);
}
#endif
