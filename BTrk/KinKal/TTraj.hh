#ifndef KinKal_TTraj_hh
#define KinKal_TTraj_hh
//
//  base class for a 1-dimensional directional path in space with time as parametric variable
//  used as part of the kinematic kalman fit
//
#include "BTrk/KinKal/Types.hh"
#include "BTrk/KinKal/TRange.hh"
namespace KinKal {

  class TTraj {
    public:
      // geometric accessors
      virtual void position(Vec4& pos) const =0; // position as a function of the input 4 vector time.
      virtual void position(double time, Vec3& pos) const =0;
      virtual void velocity(double time, Vec3& vel) const =0; // velocity vector
      virtual double speed(double time) const {
	Vec3 vel;
	velocity(time,vel);
	return sqrt(vel.Mag2());
      }
      virtual void direction(double time, Vec3& dir) const =0; // unit vector in the direction of positive time
      TTraj(TRange const& trange) : trange_(trange){}
      TTraj() {}
      TRange const& range() const { return trange_; }
      TRange& range() { return trange_; }
      bool inRange(double t) const { return trange_.inRange(t); }
    protected:
      TRange trange_;
  };

//  double TTraj::speed(double time) const {

}
#endif

