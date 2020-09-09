#ifndef KinKal_Line_hh
#define KinKal_Line_hh
//
//  Linear time-based trajectory with a constant velocity.
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/Vectors.hh"
#include "KinKal/TimeRange.hh"
namespace KinKal {
  class Line {
    public:
      // construct from a spacepoint and propagation velocity (mm/ns)
      // by default, the line has infinite unforced range
      Line(VEC4 const& p0, VEC3 const& svel, TimeRange const& range=TimeRange(),bool forcerange=false);
      Line(VEC3 const& p0, VEC3 const& svel, double tmeas, TimeRange const& range=TimeRange(),bool forcerange=false);
      // accessors
      double t0() const { return t0_; }
      VEC3 const& pos0() const { return pos0_; }
      double speed() const { return speed_; }
      // are we forcing the range?
      bool forceRange() const { return forcerange_; }
      // TOCA for a given point
      double TOCA(VEC3 point) const;

      // geometric accessors
      void position(VEC4& pos) const;
      VEC3 position(double time) const;
      VEC3 velocity(double time) const;
      VEC3 const& direction(double time) const { return dir_; }
      VEC3 const& dir() const { return dir_; }
      double speed(double time) const;
      void print(std::ostream& ost, int detail) const;
      TimeRange const& range() const { return trange_; }
      TimeRange& range() { return trange_; }
      void setRange(TimeRange const& trange) { trange_ = trange; }
      bool inRange(double time) const { return trange_.inRange(time); }

    private:
      TimeRange trange_;
      double t0_; // intial time (at pos0)
      double speed_; // signed linear velocity, translates time to distance along the trajectory (mm/nsec)
      VEC3 pos0_, dir_; // caches
      bool forcerange_; // if set, strictly enforce the range
  };
  std::ostream& operator <<(std::ostream& ost, Line const& tline);
}
#endif
