#include "KinKal/Line.hh"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  Line::Line(VEC4 const& pos0, VEC3 const& svel, TimeRange const& range,bool forcerange) : Line(pos0.Vect(), svel, pos0.T(), range, forcerange) {}
  Line::Line(VEC3 const& pos0, VEC3 const& svel, double tmeas, TimeRange const& range, bool forcerange)  : trange_(range), 
  t0_(tmeas), speed_(sqrt(svel.Mag2())), pos0_(pos0), dir_(svel.Unit()), forcerange_(forcerange) {
  }

  void Line::position(VEC4& pos) const {
    VEC3 pos3 = position(pos.T());
    pos.SetXYZT(pos3.X(),pos3.Y(),pos3.Z(),pos.T());
  }

  VEC3 Line::position(double time) const {
    if(forceRange()) range().forceRange(time);
    return pos0() + ((time-t0())*speed())*dir_;
  }

  VEC3 Line::velocity(double time) const {
    return dir_*speed();
  }

  double Line::speed(double time) const {
    return speed();
  }

  double Line::TOCA(VEC3 point) const {
    double s = (point - pos0()).Dot(dir_);
    return s/speed_ - t0();
  }

  void Line::print(std::ostream& ost, int detail) const {
    ost << " Line " <<  range()
    << " t0 " << t0_
    << " speed " << speed_
    << " intial position " << pos0_
    << " direction " << dir_ << endl;
  }

  std::ostream& operator <<(std::ostream& ost, Line const& tline) {
    tline.print(ost,0);
    return ost;
  }

}
