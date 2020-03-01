#include "KinKal/TLine.hh"
#include <iostream>
#include <stdexcept>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  vector<string> TLine::paramTitles_ = {
    "Transverse DOCA to Z Axis",
    "Azimuth of POCA"
    "Z at POCA",
    "Cos Theta",
    "Time at POCA"}; 
  vector<string> TLine::paramNames_ = {
  "D0","Phi0","Z0","CTheta","Time0"};
  std::vector<std::string> const& TLine::paramNames() { return paramNames_; }
  std::vector<std::string> const& TLine::paramTitles() { return paramTitles_; }
  std::string const& TLine::paramName(paramIndex index) { return paramNames_[static_cast<size_t>(index)];}
  std::string const& TLine::paramTitle(paramIndex index) { return paramTitles_[static_cast<size_t>(index)];}

  TLine::TLine(Vec4 const& pos0, Vec3 const& svel, TRange const& range,bool forcerange) : TLine(pos0.Vect(), svel, pos0.T(), range, forcerange) {}
  TLine::TLine(Vec3 const& pos0, Vec3 const& svel, double tmeas, TRange const& range, bool forcerange)  : TTraj(range), 
  speed_(sqrt(svel.Mag2())), pos0_(pos0), dir_(svel.Unit()), forcerange_(forcerange) {
    static const Vec3 zdir(0.0,0.0,1.0);
    double zddot = zdir.Dot(dir_);
    double stheta2 = (1.0 -zddot*zddot);
    param(cost_) = zddot;
    // separate into cases: parallel to z axis is special
    if(zddot > 1.0e-5){
      // find the POCA with the z axis; this defines the reference point
      double pddot = pos0_.Dot(dir_);
      double slen = (pos0_.Z()*zddot - pddot)/stheta2;
      auto poca = pos0_ + dir_*slen;
      param(d0_) = poca.Rho();
      param(phi0_) = atan2(poca.Y(),poca.X());
      param(z0_) = poca.Z();
      // check
      if(fabs(poca.Z()+(pddot*zddot - pos0_.Z())/stheta2) > 1e-5)
	throw std::range_error("POCA calculation failed!");
      // move the time to POCA
      param(t0_) = tmeas - slen/speed_;
    } else {
      // define parameters using the reference point
      param(d0_) = pos0_.Rho();
      param(phi0_) = atan2(pos0_.Y(),pos0_.X());
      param(z0_) = pos0_.Z();
      param(t0_) = tmeas;
    }
  }

  void TLine::position(Vec4& pos) const {
    Vec3 pos3;
    position(pos.T(),pos3);
    pos.SetXYZT(pos3.X(),pos3.Y(),pos3.Z(),pos.T());
  }

  void TLine::position(double time, Vec3& pos) const {
    if(forceRange()) range().forceRange(time);
    pos = pos0() + ((time-t0())*speed())*dir();
  }

  void TLine::velocity(double time, Vec3& vel) const {
    vel = dir()*speed();
  }

  void TLine::direction(double time, Vec3& dirvec) const {
    dirvec = dir();
  }

  double TLine::speed(double time) const {
    return speed();
  }

  double TLine::TOCA(Vec3 point) const {
    double s = (point - pos0()).Dot(dir());
    return s/speed_ - t0();
  }

  std::ostream& operator <<(std::ostream& ost, TLine const& tline) {
    ost << " TLine parameters: ";
    for(size_t ipar=0;ipar < TLine::npars_;ipar++){
      ost << TLine::paramName(static_cast<TLine::paramIndex>(ipar) ) << " : " << tline.param(ipar);
      if(ipar < TLine::npars_-1) ost << " , ";
    }
    ost <<  tline.range();
    return ost;
  }

}
