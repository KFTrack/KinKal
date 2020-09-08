#ifndef KinKal_TLine_hh
#define KinKal_TLine_hh
//
//  Linear time-based trajectory with a constant velocity.
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/Vectors.hh"
#include "KinKal/TRange.hh"
namespace KinKal {
  class TLine {
    public:
      // construct from a spacepoint and propagation velocity (mm/ns)
      // by default, the line has infinite unforced range
      TLine(Vec4 const& p0, Vec3 const& svel, TRange const& range=TRange(),bool forcerange=false);
      TLine(Vec3 const& p0, Vec3 const& svel, double tmeas, TRange const& range=TRange(),bool forcerange=false);
      // accessors
      double t0() const { return t0_; }
      Vec3 const& pos0() const { return pos0_; }
      double speed() const { return speed_; }
      // are we forcing the range?
      bool forceRange() const { return forcerange_; }
      // TOCA for a given point
      double TOCA(Vec3 point) const;

      // geometric accessors
      void position(Vec4& pos) const;
      Vec3 position(double time) const;
      Vec3 velocity(double time) const;
      Vec3 const& direction(double time) const { return dir_; }
      Vec3 const& dir() const { return dir_; }
      double speed(double time) const;
      void print(std::ostream& ost, int detail) const;
      TRange const& range() const { return trange_; }
      TRange& range() { return trange_; }
      virtual void setRange(TRange const& trange) { trange_ = trange; }
      bool inRange(double time) const { return trange_.inRange(time); }

    private:
      TRange trange_;
      double t0_; // intial time (at pos0)
      double speed_; // signed linear velocity, translates time to distance along the trajectory (mm/nsec)
      Vec3 pos0_, dir_; // caches
      bool forcerange_; // if set, strictly enforce the range
  };
  std::ostream& operator <<(std::ostream& ost, TLine const& tline);
}
#endif
