#ifndef KinKal_SensorLine_hh
#define KinKal_SensorLine_hh

//  Used as part of the kinematic Kalman fit
//
#include "KinKal/General/Vectors.hh"

namespace KinKal {
  class SensorLine {
    public:
      // construct from a spacetime point (typically the measurement position and time) and propagation velocity (mm/ns).
      SensorLine(VEC4 const& p0, VEC3 const& svel, double length);
      SensorLine(VEC3 const& p0, VEC3 const& svel, double length);
      // construct from 2 points. P0 is the measurement (near) end, p1 the far end.  Signals propagate from far to near
      SensorLine(VEC3 const& p0, VEC3 const& p1);
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
      VEC3 pos0_, dir_; // position and direction
      double length_; // line length
  };
  std::ostream& operator <<(std::ostream& ost, SensorLine const& tSensorLine);
}
#endif