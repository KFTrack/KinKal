#include "KinKal/TLine.hh"
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  TLine::TLine(Vec4 const& pos0, Vec3 const& svel, TRange const& range,bool forcerange) : TLine(pos0.Vect(), svel, pos0.T(), range, forcerange) {}
  TLine::TLine(Vec3 const& pos0, Vec3 const& svel, double tmeas, TRange const& range, bool forcerange)  : trange_(range), 
  t0_(tmeas), speed_(sqrt(svel.Mag2())), pos0_(pos0), dir_(svel.Unit()), forcerange_(forcerange) {
  }

  void TLine::position(Vec4& pos) const {
    Vec3 pos3 = position(pos.T());
    pos.SetXYZT(pos3.X(),pos3.Y(),pos3.Z(),pos.T());
  }

  Vec3 TLine::position(double time) const {
    if(forceRange()) range().forceRange(time);
    return pos0() + ((time-t0())*speed())*dir_;
  }

  Vec3 TLine::velocity(double time) const {
    return dir_*speed();
  }

  double TLine::speed(double time) const {
    return speed();
  }

  double TLine::TOCA(Vec3 point) const {
    double s = (point - pos0()).Dot(dir_);
    return s/speed_ - t0();
  }

  void TLine::print(std::ostream& ost, int detail) const {
    ost << " TLine " <<  range()
    << " t0 " << t0_
    << " speed " << speed_
    << " intial position " << pos0_
    << " direction " << dir_ << endl;
  }

  std::ostream& operator <<(std::ostream& ost, TLine const& tline) {
    tline.print(ost,0);
    return ost;
  }

}
