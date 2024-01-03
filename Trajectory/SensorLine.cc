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

  SensorLine::SensorLine(VEC4 const& mpos, VEC3 const& svel, double length ) : SensorLine(mpos.Vect(), mpos.T(), svel, length) {}

// note the measurement position is at the END of the directed line segment, since that points in the direction of signal propagation
  SensorLine::SensorLine(VEC3 const& mpos, double mtime , VEC3 const& svel, double length )  :
    mtime_(mtime), d2t_(new ConstantDistanceToTime(svel.R())), lineseg_(mpos-svel.unit()*length, mpos) {}

  SensorLine::SensorLine(VEC3 const& mpos, VEC3 const& endpos, double mtime, double speed ) : mtime_(mtime),
  d2t_(new ConstantDistanceToTime(speed)), lineseg_(endpos,mpos) {}

  SensorLine::SensorLine(VEC3 const& mpos, double mtime, VEC3 const& svel, double length, std::shared_ptr<DistanceToTime> d) : mtime_(mtime),
  d2t_(d), lineseg_(mpos-svel.unit()*length, mpos) {}

  VEC3 SensorLine::position3(double time) const {
    return lineseg_.position(d2t_->distance(time - mtime_));
  }

  VEC4 SensorLine::position4(double time) const {
    VEC3 pos3 = position3(time);
    return VEC4(pos3.X(),pos3.Y(),pos3.Z(),time);
  }

  VEC3 SensorLine::velocity(double time) const {
    return direction(time)*speed(time);
  }

  double SensorLine::timeTo(VEC3 const& point) const {
    double s = lineseg_.distanceTo(point);
    return s/speed(s) - measurementTime();
  }

  void SensorLine::print(std::ostream& ost, int detail) const {
    ost << *this << endl;
  }

  std::ostream& operator <<(std::ostream& ost, SensorLine const& sline) {
    ost << "SensorLine " << sline.line() << " measurement time " << sline.measurementTime();
    return ost;
  }

}
