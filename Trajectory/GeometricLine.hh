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
      // construct from 2 points. P0 is the signal origin position (start position), p1 the measurement position(end position).  Signals propagate from start to end.
      GeometricLine(VEC3 const& p0, VEC3 const& p1);
      // accessors
      // signal ends at pos0
      VEC3 startPosition() const { return pos0_ - length_*dir_; }
      VEC3 const& endPosition() const { return pos0_ ; }
      double length() const { return length_; }
      VEC3 const& direction() const { return dir_; }
      // Distance of Closest Approach to a point
      double DOCA(VEC3 const& point) const;
      // geometric accessors
      VEC3 position3(double distance) const;
      void print(std::ostream& ost, int detail) const;

    private:
      VEC3 pos0_, dir_; // signal physical origin position and signal propagation direction
      double length_; // distance from signal origin to measurement (electronics) point
  };
  std::ostream& operator <<(std::ostream& ost, GeometricLine const& tGeometricLine);
}
#endif
