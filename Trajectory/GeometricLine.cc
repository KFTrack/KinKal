#include "KinKal/Trajectory/GeometricLine.hh"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  GeometricLine::GeometricLine(VEC4 const& p0, VEC3 const& svel, double length) : GeometricLine(p0.Vect(), svel, length) {}
  GeometricLine::GeometricLine(VEC3 const& p0, VEC3 const& svel, double length)  : 
    pos0_(p0), dir_(svel.Unit()), length_(length) {}
  
  // flipped
  // GeometricLine::GeometricLine(VEC3 const& p0, VEC3 const& p1) : pos0_(p0), dir_((p0-p1).Unit()), length_((p1-p0).R()) {}
  GeometricLine::GeometricLine(VEC3 const& p0, VEC3 const& p1) : pos0_(p0), dir_((p1-p0).Unit()), length_((p1-p0).R()) {}

  VEC3 GeometricLine::position3(double distance) const {
    return pos0_ + distance*dir_;
  }

  double GeometricLine::DOCA(VEC3 const& point) const {
    // flipped
    // return (point - pos0_).Dot(dir_);
    return (pos0_ - point).Dot(dir_);
  }

  void GeometricLine::print(std::ostream& ost, int detail) const {
    ost << " GeometricLine, initial position " << pos0_
    << " direction " << dir_
    << " length " << length_ << endl;
  }

  std::ostream& operator <<(std::ostream& ost, GeometricLine const& tGeometricLine) {
    tGeometricLine.print(ost,0);
    return ost;
  }

}
