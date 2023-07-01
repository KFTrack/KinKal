#ifndef KinKal_GeometricLine_hh
#define KinKal_GeometricLine_hh

//  Used as part of the kinematic Kalman fit
//
#include "KinKal/General/Vectors.hh"

namespace KinKal {
  class GeometricLine {
    public:
      // construct from a spacetime point (typically the measurement position and time) and propagation velocity (mm/ns).
      GeometricLine(VEC4 const& p0, VEC3 const& svel, double length);
      GeometricLine(VEC3 const& p0, VEC3 const& svel, double length);
      // construct from 2 points. P0 is the measurement (near) end, p1 the far end.  Signals propagate from far to near
      GeometricLine(VEC3 const& p0, VEC3 const& p1);
      // accessors
      // signal ends at pos0
      // flipped
      // VEC3 startPosition() const { return pos0_ - length_*dir_; }
      // VEC3 const& endPosition() const { return pos0_; }
      VEC3 const& startPosition() const { return pos0_; }
      VEC3 endPosition() const { return pos0_ + length_*dir_; }
      double length() const { return length_; }
      VEC3 const& direction() const { return dir_; }
      // Distance of Closest Approach to a point
      double DOCA(VEC3 const& point) const;
      // geometric accessors
      VEC3 position3(double distance) const;
      void print(std::ostream& ost, int detail) const;

    private:
      VEC3 pos0_, dir_; // position and direction
      double length_; // line length
  };
  std::ostream& operator <<(std::ostream& ost, GeometricLine const& tGeometricLine);
}
#endif