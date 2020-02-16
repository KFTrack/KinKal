#include "BTrk/KinKal/TLine.hh"
#include "BTrk/KinKal/Constants.hh"
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

  TLine::TLine(Vec3 const& p0, Vec3 const& svel, double tmeas, TRange const& range) : TTraj(range) {
    speed_ = sqrt(svel.Mag2());
    static const Vec3 zdir(0.0,0.0,1.0);
    auto sdir = svel.Unit();
    double zsdot = zdir.Dot(sdir);
    double stheta2 = (1.0 -zsdot*zsdot);
    param(cost_) = zsdot;
    // separate into cases: parallel to z axis is special
    if(zsdot > 1.0e-5){
      // find the POCA with the z axis; this defines the reference point
      double psdot = p0.Dot(sdir);
      double slen = (p0.Z()*zsdot - psdot)/stheta2;
      auto poca = p0 + sdir*slen;
      param(d0_) = poca.Rho();
      param(phi0_) = atan2(poca.Y(),poca.X());
      param(z0_) = poca.Z();
    // check
      if(fabs(poca.Z()+(psdot*zsdot - p0.Z())/stheta2) > 1e-5)
	throw std::range_error("POCA calculation failed!");
    // move the time to POCA
      param(t0_) = tmeas - slen/speed_;
    } else {
    // define parameters using the reference point
      param(d0_) = p0.Rho();
      param(phi0_) = atan2(p0.Y(),p0.X());
      param(z0_) = p0.Z();
      param(t0_) = tmeas;
    }
    // set the cache values
    double sint = sinTheta();
    dir_.SetXYZ(sint*sin(phi0()), -sint*cos(phi0()), cost());
    
  }

  void TLine::position(Vec4& pos) const {
    Vec3 pos3;
    position(pos.T(),pos3);
    pos.SetXYZT(pos3.X(),pos3.Y(),pos3.Z(),pos.T());
  }

  void TLine::position(double time, Vec3& pos) const {
    if(forceRange()) range().forceRange(time);
    pos = p0() + (time-t0())*dir()*speed()
  }

  void TLine::velocity(double time, Vec3& vel) const {
    vel = dir()*speed();
  }

  void TLine::direction(double time, Vec3& dir) const {
    dir = dir();
  }

  double TLine::speed(double time) const {
    return speed();
  }

  double TLine::TOCA(Vec3 point) const {
    double s = (point - p0()).Dot(dir());
    return s/speed_ - t0();
  }

}
