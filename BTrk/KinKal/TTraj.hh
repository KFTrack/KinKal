#ifndef kinkal_ttraj_hh
#define kinkal_ttraj_hh
//
//  base class for a trajectory (1-dimensional path in space with time as parametric variable
//  used as part of the kinematic kalman fit
//
#include "BTrk/KinKal/Types.hh"
namespace KinKal {
  class TTraj {
    public:
      // geometric accessors
      virtual void position(Vec4& pos) const =0; // position as a function of the input 4 vector time.
      virtual void position(double time, Vec3& pos) const =0;
      virtual void velocity(double time, Vec3& vel) const =0; // velocity vector in mm/ns
      virtual void direction(double time, Vec3& dir) const =0; // unit vector in the direction of positive time
      virtual ~TTraj() = 0;
  };
}
#endif

