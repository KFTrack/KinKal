#include "KinKal/Trajectory/SensorLine.hh"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  SensorLine::SensorLine(VEC4 const& p0, VEC3 const& svel, double length) : SensorLine(p0.Vect(), svel, length) {}
  SensorLine::SensorLine(VEC3 const& p0, VEC3 const& svel, double length )  : 
    pos0_(p0), dir_(svel.Unit()), length_(length) {}
  SensorLine::SensorLine(VEC3 const& p0, VEC3 const& p1) : pos0_(p0), dir_((p0-p1).Unit()), length_((p1-p0).R()) {}

  VEC3 SensorLine::position3(double distance) const {
    return pos0_ + distance*dir_;
  }

  double SensorLine::DOCA(VEC3 const& point) const {
    return (point - pos0_).Dot(dir_);
  }

  void SensorLine::print(std::ostream& ost, int detail) const {
    ost << " SensorLine, initial position " << pos0_
    << " direction " << dir_
    << " length " << length_ << endl;
  }

  std::ostream& operator <<(std::ostream& ost, SensorLine const& tSensorLine) {
    tSensorLine.print(ost,0);
    return ost;
  }

}
