#include "KinKal/General/CylBFieldMap.hh"

namespace KinKal{
  using VEC = std::vector<double>;
  using MAT = std::vector<std::vector<double>>;
  
  CylBFieldMap::CylBFieldMap(std::string const& file)
    : B_dat(file), Br_interp(B_dat.r_, B_dat.z_, B_dat.Br_),
    Bz_interp(B_dat.r_, B_dat.z_, B_dat.Bz_) {}

  VEC3 CylBFieldMap::fieldVect(VEC3 const& position) const {
    // calcualte interpolated field value at position, in cartesian coordinates
    double Bz = Bz_interp.interp(position.Rho(), position.Z());
    double Br = Br_interp.interp(position.Rho(), position.Z());
    // handle case where R=0 -- R=X
    if (position.Rho() < 1e-6) {
      return VEC3(Br, 0., Bz);
    }
    else {
      return VEC3(Br*position.X()/position.Rho(), Br*position.Y()/position.Rho(), Bz);
    }
  }

  CylBFieldMap::Grad CylBFieldMap::fieldGrad(VEC3 const& position) const {
    // calculate gradient matrix at position, in cartesian coordinates
    // calculate needed field and gradient values
    double Br = Br_interp.interp(position.Rho(), position.Z());
    VEC dBz = Bz_interp.gradient(position.Rho(), position.Z());
    VEC dBr = Br_interp.gradient(position.Rho(), position.Z());
    // vector to fill, will instantiate final Grad
    VEC dBi(9);
    // handle case where R=0 -- R=X, dR/dX=1, dR/dY=0
    if (position.Rho() < 1e-6) {
      dBi = VEC{dBr[0], 0., dBr[1],
      0.,0.,0.,
      dBz[0], 0., dBz[1]};
    }
    else {
      // store values that are calculated multiple times
      double xr=position.X()/position.Rho(), yr=position.Y()/position.Rho();
      double xryr = xr*yr, xrxr=xr*xr, yryr=yr*yr;
      double rinv = 1/position.Rho();
      dBi = VEC{dBr[0]*xrxr + Br * rinv*(1 - xrxr),
        dBr[0]*xryr + Br * (- xryr * rinv),
        xr * dBr[1],
        dBr[0]*xryr + Br * (- xryr * rinv),
        dBr[0]*yryr + Br * rinv*(1 - yryr),
        yr * dBr[1],
        dBz[0]*xr,dBz[0]*yr,dBz[1]};
    }
    Grad retval(dBi.begin(), dBi.end());
    return retval;
  }

  VEC3 CylBFieldMap::fieldDeriv(VEC3 const& position, VEC3 const& velocity) const {
    Grad grad = fieldGrad(position);
    VEC3 retval(grad[0][0]*velocity.X()+grad[0][1]*velocity.Y()+grad[0][2]*velocity.Z(),
      grad[1][0]*velocity.X()+grad[1][1]*velocity.Y()+grad[1][2]*velocity.Z(),
      grad[2][0]*velocity.X()+grad[2][1]*velocity.Y()+grad[2][2]*velocity.Z());
    return retval;
  }

  bool CylBFieldMap::inRange(VEC3 const& position) const {
    // check if position is in range of the map values
    bool rrange = (position.R() >= 0.) && (position.R() <= rMax());
    bool zrange = (position.Z() >= zMin()) && (position.Z() <= zMax());
    return rrange && zrange;
  }
}
