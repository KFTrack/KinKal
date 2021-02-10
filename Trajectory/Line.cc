#include "KinKal/Trajectory/Line.hh"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  Line::Line(VEC4 const& pos0, VEC3 const& svel, double length ) : Line(pos0.Vect(), pos0.T(), svel, length) {}
  Line::Line(VEC3 const& pos0, double tmeas , VEC3 const& svel, double length )  : 
    pos0_(pos0), dir_(svel.Unit()), t0_(tmeas), speed_(sqrt(svel.Mag2())), length_(length) {}
  Line::Line(VEC3 const& p0, VEC3 const& p1, double t0, double speed ) : pos0_(p0), dir_((p0-p1).Unit()), t0_(t0), speed_(speed), length_((p1-p0).R()) {}

  VEC3 Line::position3(double time) const {
    return pos0_ + ((time-t0_)*speed_)*dir_;
  }

  VEC4 Line::position4(double time) const {
    VEC3 pos3 = position3(time);
    return VEC4(pos3.X(),pos3.Y(),pos3.Z(),time);
  }

  VEC3 Line::velocity(double time) const {
    return dir_*speed_;
  }

  double Line::TOCA(VEC3 point) const {
    double s = (point - pos0_).Dot(dir_);
    return s/speed_ - t0_;
  }

  void Line::print(std::ostream& ost, int detail) const {
    ost << " Line, intial position " << pos0_
    << " t0 " << t0_
    << " direction " << dir_
    << " speed " << speed_ << endl;
  }

  std::ostream& operator <<(std::ostream& ost, Line const& tline) {
    tline.print(ost,0);
    return ost;
  }

}
