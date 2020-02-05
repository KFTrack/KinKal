#ifndef KinKal_TTraj_hh
#define KinKal_TTraj_hh
//
//  base class for a 1-dimensional directional path in space with time as parametric variable
//  used as part of the kinematic kalman fit
//
#include "BTrk/KinKal/Types.hh"
namespace KinKal {

// simple struct to describe a time range
  struct TRange {
    std::array<double,2> range_; // range of times
    TRange() : range_{1.0,-1.0} {} // initialize to have infinite range
    bool inRange(double t) const { return (range_[0] > range_[1]) ||
      (t > range_[0] && t < range_[1]); }
  };

  class TTraj {
    public:
      // geometric accessors
      virtual void position(Vec4& pos) const =0; // position as a function of the input 4 vector time.
      virtual void position(double time, Vec3& pos) const =0;
      virtual void velocity(double time, Vec3& vel) const =0; // velocity vector
      virtual void direction(double time, Vec3& dir) const =0; // unit vector in the direction of positive time
      TTraj(TRange const& trange) : trange_(trange){}
      TTraj() {}
      bool inRange(double t) const { return trange_.inRange(t); }
    private:
      TRange trange_;
  };
}
#endif

