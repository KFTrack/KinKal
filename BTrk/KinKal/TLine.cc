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
    pars_.vec()[cost_] = zsdot;
    // separate into cases: parallel to z axis is special
    if(zsdot > 1.0e-5){
      // find the POCA with the z axis; this defines the reference point
      double psdot = p0.Dot(sdir);
      double slen = (p0.Z()*zsdot - psdot)/stheta2;
      auto poca = p0 + sdir*slen;
      pars_.vec()[d0_] = poca.Rho();
      pars_.vec()[phi0_] = atan2(poca.Y(),poca.X());
      pars_.vec()[z0_] = poca.Z();
    // check
      if(fabs(poca.Z()+(psdot*zsdot - p0.Z())/stheta2) > 1e-5)
	throw std::range_error("POCA calculation failed!");
    // move the time to POCA
      pars_.vec()[t0_] = tmeas - slen/speed_;
    } else {
    // define parameters using the reference point
      pars_.vec()[d0_] = p0.Rho();
      pars_.vec()[phi0_] = atan2(p0.Y(),p0.X());
      pars_.vec()[z0_] = p0.Z();
      pars_.vec()[t0_] = tmeas;
    }
  }

  void TLine::position(Vec4& pos) const {
    Vec3 pos3;
    position(pos.T(),pos3);
    pos.SetXYZT(pos3.X(),pos3.Y(),pos3.Z(),pos.T());
  }

  void TLine::pos0( Vec3& pos) const {
    pos = Vec3(d0()*cos(phi0()), d0()*sin(phi0()), z0());
  }

  void TLine::position(double time, Vec3& pos) const {
    Vec3 p0; pos0(p0);
    Vec3 vel; velocity(time,vel);
    pos = p0 + (time-t0())*vel;
  }


  void TLine::velocity(double time, Vec3& vel) const {
    direction(time,vel);
    vel *= speed_;
  }

  void TLine::direction(double time, Vec3& dir) const {
    double sint = sinTheta();
    dir.SetXYZ(sint*sin(phi0()), -sint*cos(phi0()), cost());
  }

  double TLine::speed(double time) const {
    return speed_;
  }

}
