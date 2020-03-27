#ifndef KinKal_TTraj_hh
#define KinKal_TTraj_hh
//
//  base class for a 1-dimensional directional path in space with time as parametric variable
//  used as part of the kinematic kalman fit
//
#include "KinKal/Vectors.hh"
#include "KinKal/TRange.hh"
namespace KinKal {
  class BField;
  class TTraj {
    public:
      // geometric accessors
      virtual void position(Vec4& pos) const =0; // position as a function of the input 4 vector time.
      virtual void position(double time, Vec3& pos) const =0;
      virtual void velocity(double time, Vec3& vel) const =0; // velocity vector
      virtual double speed(double time) const =0;
      virtual void direction(double time, Vec3& dir) const =0; // unit vector in the direction of positive time
      // constructors, etc
      TTraj(TRange const& trange) : trange_(trange){}
      TTraj() {}
      virtual ~TTraj() {}
      // range is an intrinsic part
      TRange const& range() const { return trange_; }
      TRange& range() { return trange_; }
      virtual void setRange(TRange const& trange) { trange_ = trange; }
      bool inRange(double t) const { return trange_.inRange(t); }
    protected:
      TRange trange_;
  };
}
#endif

