#ifndef KinKal_Line_hh
#define KinKal_Line_hh
//
//  Linear time-based trajectory with a constant velocity.
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/General/Vectors.hh"
#include "KinKal/General/TimeRange.hh"
namespace KinKal {
  class Line {
    public:
      // construct from a spacetime point (typically the measurement position and time) and propagation velocity (mm/ns).
      Line(VEC4 const& p0, VEC3 const& svel, double length);
      Line(VEC3 const& p0, double t0, VEC3 const& svel, double length);
      // construct from 2 points plus timing information.  P0 is the measurement (near) end, p1 the far end.  Signals propagate from far to near
      Line(VEC3 const& p0, VEC3 const& p1, double t0, double speed );
      // accessors
      double t0() const { return t0_; }
      double& t0() { return t0_; } // detector updates need to refine t0
      VEC3 const& startPosition() const { return pos0_; }
      VEC3 endPosition() const { return pos0_ + length_*dir_; }
      double speed() const { return speed_; }
      double speed(double time) const { return speed_; }
      double length() const { return length_; }
      VEC3 const& direction() const { return dir_; }
      // TOCA to a point
      double TOCA(VEC3 const& point) const;
      // geometric accessors
      VEC3 position3(double time) const;
      VEC4 position4(double time) const;
      VEC3 velocity(double time) const;
      VEC3 const& direction(double time) const { return dir_; }
      void print(std::ostream& ost, int detail) const;
      TimeRange range() const { return TimeRange(t0_,t0_ +length_/speed_); }

    private:
      VEC3 pos0_, dir_; // position and direction
      double t0_; // intial time (at pos0)
      double speed_; // signed linear velocity, translates time to distance along the trajectory (mm/nsec)
      double length_; // line length
  };
  std::ostream& operator <<(std::ostream& ost, Line const& tline);
}
#endif
