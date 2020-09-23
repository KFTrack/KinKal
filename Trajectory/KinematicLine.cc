/*
  KinematicLine is the Linear Trajectory Specialization of KTRAJ - the kinematic trajectory.
  Original Author: S Middleton 2020
*/

#include "KinKal/Trajectory/KinematicLine.hh"
#include "KinKal/Detector/BFieldMap.hh"
#include "KinKal/Trajectory/POCAUtil.hh"
#include "Math/AxisAngle.h"
#include "Math/VectorUtil.h"
#include <math.h>
#include <stdexcept>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  typedef ROOT::Math::SVector<double,3> SVEC3;
  const vector<string> KinematicLine::paramTitles_ = {
      "Transverse DOCA to Z Axis (d_{0})", "Azimuth of POCA (#phi_{0})",
      "Z at POCA (z_{0})", "Cos #theta", "Time at POCA (t_{0})", "Momentum"};

  const vector<string> KinematicLine::paramNames_ = {"d_{0}", "#phi_{0}", "z_{0}",
                                       "cos(#theta)", "t_{0}", "mom"};

  const vector<string> KinematicLine::paramUnits_ = {"mm", "radians", "mm", "", "ns","MeV/c"};

  std::vector<std::string> const &KinematicLine::paramUnits() { return paramUnits_; }
  std::vector<std::string> const &KinematicLine::paramNames() { return paramNames_; }
  std::vector<std::string> const &KinematicLine::paramTitles() { return paramTitles_; }

  std::string const &KinematicLine::paramName(ParamIndex index) {
    return paramNames_[static_cast<size_t>(index)];
  }
  std::string const &KinematicLine::paramTitle(ParamIndex index) {
    return paramTitles_[static_cast<size_t>(index)];
  }
  std::string const &KinematicLine::paramUnit(ParamIndex index) {
    return paramUnits_[static_cast<size_t>(index)];
  }

  const string KinematicLine::trajName_("KinematicLine");
  string const &KinematicLine::trajName() { return trajName_; }

 KinematicLine::KinematicLine( VEC4 const& pos0, MOM4 const& mom0, int charge, double bnom, TimeRange const& range) : KinematicLine(pos0,mom0,charge,VEC3(0.0,0.0,bnom),range) {}

  KinematicLine::KinematicLine(VEC4 const &pos0, MOM4 const &mom0, int charge, VEC3 const &bnom,
  TimeRange const &trange) :  bnom_(bnom), mass_(mom0.M()), charge_(charge), trange_(trange)
  {
    double mommag = mom0.R();
    double speed = ( mommag/ mom0.E()) * CLHEP::c_light;
    VEC3 dir = mom0.Vect().Unit();
    
    static const VEC3 zdir(0.0, 0.0, 1.0);
    static const VEC3 zpos(0.0, 0.0, 0.0);

// calculate POCA to the Z axis.  This is the reference for the parameters
    POCAUtil poca(pos0.Vect(), dir, zpos, zdir);
    double flen = dir.Dot(pos0.Vect()-poca.point1()); // flight length from reference to POCA
    VEC3 pca = poca.point1()-poca.point2(); // vector from Z axis to POCA
// doca is signed by the angular momentum around the Z axis
    double amsign = copysign(1.0, -(zdir.Cross(pca)).Dot(dir));
    param(d0_) = amsign*poca.dca(); // dca2d and dca are the same for POCA to the Z axis 
    param(phi0_) = dir.Phi(); // same as at POCA
    param(z0_) = poca.point1().Z();
    param(cost_) = dir.Z();
    param(t0_) = pos0.T() - flen/speed;
    param(mom_) = mommag;
  }

  /*
  KinematicLine can take in Momentum externally as a 4-vector or calculate it based. You
  can initialize the line with an origin (pos0) or the trajectory parameters
  (pdata)
  */

  KinematicLine::KinematicLine(Parameters const &pdata, KinematicLine const &other) : KinematicLine(other) {
    pars_ = pdata;
  }

  KinematicLine::KinematicLine(ParticleState const& pstate, int charge, VEC3 const& bnom, TimeRange const& range) :
    KinematicLine(pstate.position4(),pstate.momentum4(), charge,bnom,range) 
  {}

  KinematicLine::KinematicLine(ParticleStateMeasurement const& pstate, int charge, VEC3 const& bnom, TimeRange const& range) :
  KinematicLine(pstate.stateVector(),charge,bnom,range) {
  // derive the parameter space covariance from the global state space covariance
    PSMAT dpds = dPardState(pstate.stateVector().time());
    pars_.covariance() = ROOT::Math::Similarity(dpds,pstate.stateCovariance());
  }

  KinematicLine::KinematicLine(KinematicLine const& other, VEC3 const& bnom, double trot) : KinematicLine(other) {
  }

  VEC3 KinematicLine::position3(double time) const {
    return (pos0() + flightLength(time) * dir());
  }

  VEC4 KinematicLine::position4(double time) const {
    VEC3 temp = position3(time);
    return VEC4(temp.X(), temp.Y(), temp.Z(), time);
  }

  VEC3 KinematicLine::momentum3(double tval) const {
    return direction(tval)*mom();
  }

  MOM4 KinematicLine::momentum4(double tval) const {
    VEC3 mom3 = momentum3(tval);
    return MOM4(mom3.X(),mom3.Y(),mom3.Z(),mass_);
  }

  ParticleState KinematicLine::state(double time) const {
    return ParticleState(position4(time),momentum4(time));
  }

  ParticleStateMeasurement KinematicLine::measurementState(double time) const {
    // express the parameter space covariance in global state space
    PSMAT dsdp = dStatedPar(time);
    return ParticleStateMeasurement(state(time),ROOT::Math::Similarity(dsdp,pars_.covariance()));
  }

  /*
     The effects for changes in 2 perpendicular directions (theta1 = theta and
     theta2 = phi()*sin(theta) can sometimes be added, as scattering in these
  are uncorrelated. These axes are track specific. as cosmics are not always
  coming along the same track direction it is necessary to have difference
  parameterization than that used for the helix case.
  alt dir = a test with the "BTrk parameterization" - just changes signs due to
  swithc in cos<->sin
  */
  VEC3 KinematicLine::direction(double t, MomBasis::Direction mdir) const {

    switch (mdir) {
    case MomBasis::perpdir_: // purely polar change theta 1 = theta
      return VEC3(cosTheta() * cosPhi0(), cosTheta() * sinPhi0(), -1 * sinTheta());
      break;
    case MomBasis::phidir_: // purely transverse theta2 = -phi()*sin(theta)
      return VEC3(-sinPhi0(), cosPhi0(), 0.0);
      break;
    case MomBasis::momdir_: // along momentum
      return dir();
      break;
    default:
      throw std::invalid_argument("Invalid direction");
    }
  }

  PSMAT KinematicLine::dPardState(double time) const{
  // aggregate state from separate X and M derivatives; parameter space is row
    DPDV dPdX = dPardX(time);
    DPDV dPdM = dPardM(time);
    PSMAT dpds;
    dpds.Place_at(dPdX,0,0);
    dpds.Place_at(dPdM,0,3);
    return dpds;
  }

  PSMAT KinematicLine::dStatedPar(double time) const {
  // aggregate state from separate X and M derivatives; parameter space is column
    DVDP dXdP = dXdPar(time);
    DVDP dMdP = dMdPar(time);
    PSMAT dsdp;
    dsdp.Place_at(dXdP,0,0);
    dsdp.Place_at(dMdP,3,0);
    return dsdp;
  }

  DVDP KinematicLine::dXdPar(double time) const { 
    double deltat = time-t0();
    double sinT = sinTheta();
    double cotT = 1.0/tanTheta();
    double cosT = cosTheta();
    double sinF = sin(phi0());
    double cosF = cos(phi0());
    double spd = speed();
    double gam = gamma();
    SVEC3 dX_dd0(-sinF, cosF, 0.0);
    SVEC3 dX_dphi0 = (sinT*spd*deltat)*dX_dd0 - d0()*SVEC3(cosF,sinF,0.0);
    SVEC3 dX_dz0 (0.0,0.0,1.0);
    SVEC3 dX_dcost = (spd*deltat)*SVEC3(-cotT*cosF,-cotT*sinF,1.0);
    SVEC3 dX_dt0 = -spd*SVEC3(sinT*cosF,sinT*sinF,cosT);
    SVEC3 dX_dmom = -(deltat/(gam*gam*mom()))*dX_dt0;

    DVDP dXdP;
    dXdP.Place_in_col(dX_dd0,0,d0_);
    dXdP.Place_in_col(dX_dphi0,0,phi0_);
    dXdP.Place_in_col(dX_dz0,0,z0_);
    dXdP.Place_in_col(dX_dcost,0,cost_);
    dXdP.Place_in_col(dX_dt0,0,t0_);
    dXdP.Place_in_col(dX_dmom,0,mom_);
    return dXdP;
  }
  DVDP KinematicLine::dMdPar(double time) const {
    double sinT = sinTheta();
    double cotT = 1.0/tanTheta();
    double cosT = cosTheta();
    double sinF = sin(phi0());
    double cosF = cos(phi0());

// note: dMdd0 = dMdz0 = dMdt0 = 0
    SVEC3 dM_dphi0 = (mom()*sinT)*SVEC3(-sinF,cosF,0);
    SVEC3 dM_dcost = mom()*SVEC3(-cotT*cosF,-cotT*sinF,1.0);
    SVEC3 dM_dmom = SVEC3(sinT*cosF,sinT*sinF,cosT);

    DVDP dMdP;
    dMdP.Place_in_col(dM_dphi0,0,phi0_);
    dMdP.Place_in_col(dM_dcost,0,cost_);
    dMdP.Place_in_col(dM_dmom,0,mom_);
    return dMdP;
  }

  DPDV KinematicLine::dPardX(double time) const {
    double sinT = sinTheta();
    double cotT = 1.0/tanTheta();
    double sinF = sin(phi0());
    double cosF = cos(phi0());
    double E = energy();
// note: dCosTdX = dPhi0dX = dmom_dX = 0
    SVEC3 dd0_dX = SVEC3(-sinF,cosF,0.0);
    SVEC3 dz0_dX = SVEC3(-cotT*cosF,-cotT*sinF,1.0);
    SVEC3 dt0_dX = -E*SVEC3(cosF,sinF,0.0)/(mom()*sinT*CLHEP::c_light);
    DPDV dPdX;
    dPdX.Place_in_row(dd0_dX,d0_,0);
    dPdX.Place_in_row(dz0_dX,z0_,0);
    dPdX.Place_in_row(dt0_dX,t0_,0);

    return dPdX;
  }

  DPDV KinematicLine::dPardM(double time) const { 
    double sinT = sinTheta();
    double cosT = cosTheta();
    double cotT = 1.0/tanTheta();
    double sinF = sin(phi0());
    double cosF = cos(phi0());
    double cos2F = cosF*cosF-sinF*sinF;
    double sin2F = 2*cosF*sinF;
    VEC3 momv = momentum3(time);
    VEC3 pos = position3(time);
    static VEC3 zdir(0.0,0.0,1.0);
    VEC3 momt = VectorUtil::PerpVector(momv,zdir);
    double momt2 = momt.Mag2();
    double xmt = pos.Dot(momt);
    double E = energy();
    SVEC3 dmom_dM(sinT*cosF, sinT*sinF, cosT);
    SVEC3 dcost_dM = (1.0/mom())*(SVEC3(0.0,0.0,1.0) - cosT*dmom_dM);
    SVEC3 dphi0_dM = (1.0/(mom()*sinT))*SVEC3(-sinF,cosF,0.0);
    SVEC3 dt0_dM = (1.0/(momt2*CLHEP::c_light))*(
	xmt*( (2.0*E/momt2)*SVEC3(momv.X(),momv.Y(),0.0) 
	  - (1.0/E)*( SVEC3(momv.X(),momv.Y(),momv.Z())) )
	-  E*SVEC3(pos.X(),pos.Y(),0.0));
    SVEC3 dz0_dM = (1.0/(mom()*sinT))*(SVEC3(cotT*(cos2F*pos.X() + sin2F*pos.Y()),cotT*(sin2F*pos.X()-cos2F*pos.Y()),-cosF*pos.X()-sinF*pos.Y()));
    SVEC3 dd0_dM  = ( xmt/momt2 )* SVEC3(sinF, -cosF, 0.0);
    DPDV dPdM;
    dPdM.Place_in_row(dmom_dM,mom_,0);
    dPdM.Place_in_row(dcost_dM,cost_,0);
    dPdM.Place_in_row(dphi0_dM,phi0_,0);
    dPdM.Place_in_row(dz0_dM,z0_,0);
    dPdM.Place_in_row(dt0_dM,t0_,0);
    dPdM.Place_in_row(dd0_dM,d0_,0);
    return dPdM;
  }

  // derivatives of momentum projected along the given basis WRT the parameters
  DVEC KinematicLine::momDeriv(double time, MomBasis::Direction mdir) const {
    DPDV dPdM = dPardM(time);
    auto dir = direction(time,mdir);
    double mommag = momentum(time);
    return mommag*(dPdM*SVEC3(dir.X(), dir.Y(), dir.Z()));
  }

  void KinematicLine::print(ostream &ost, int detail) const {
    auto perr = params().covariance().Diagonal();
    ost << " KinematicLine " << range() << " parameters: ";
    for (size_t ipar = 0; ipar < KinematicLine::npars_; ipar++) {
      ost << KinematicLine::paramName(static_cast<ParamIndex>(ipar)) << " "
          << paramVal(ipar) << " +- " << perr(ipar);
      if (ipar < KinematicLine::npars_ - 1)
        ost << " ";
    }
    ost << " with rotation around Bnom " << bnom_ << endl;
  }

  ostream &operator<<(ostream &ost, KinematicLine const &lhel) {
    lhel.print(ost, 0);
    return ost;
  }

} // namespace KinKal
