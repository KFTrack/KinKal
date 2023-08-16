#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Trajectory/SensorLine.hh"
#include "KinKal/Trajectory/ConstantDistanceToTime.hh"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
#include <memory>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  Line::Line(VEC4 const& pos0, VEC3 const& svel, double length ) : Line(pos0.Vect(), pos0.T(), svel, length) {}
  Line::Line(VEC3 const& pos0, double tmeas , VEC3 const& svel, double length )  :
    t0_(tmeas), d2t_(new ConstantDistanceToTime(sqrt(svel.Mag2()))), sline_(pos0, svel, length) {}
  Line::Line(VEC3 const& p0, VEC3 const& p1, double t0, double speed ) : t0_(t0), d2t_(new ConstantDistanceToTime(speed)), sline_(p0, p1) {}
  Line::Line(VEC3 const& p0, double length, VEC3 const& svel, double t0, std::shared_ptr<DistanceToTime> d) : t0_(t0), d2t_(d), sline_(p0, svel, length) {}

  VEC3 Line::position3(double time) const {
    return sline_.position3(d2t_->distance(time - t0_));
  }

  VEC4 Line::position4(double time) const {
    VEC3 pos3 = position3(time);
    return VEC4(pos3.X(),pos3.Y(),pos3.Z(),time);
  }

  VEC3 Line::velocity(double time) const {
    return direction(time)*speed(time);
  }

  double Line::TOCA(VEC3 const& point) const {
    double s = sline_.DOCA(point);
    return s/speed(s) - t0();
  }

  void Line::print(std::ostream& ost, int detail) const {
    ost << " Line, initial position " << startPosition()
    << " t0 " << t0()
    << " direction " << direction() << endl;
  }

  std::ostream& operator <<(std::ostream& ost, Line const& tline) {
    tline.print(ost,0);
    return ost;
  }

}
