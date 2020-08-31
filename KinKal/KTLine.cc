/*
  KTLine is the Linear Trajectory Specialization of KTRAJ - the kinematic trajectory.
  Original Author: S Middleton 2020
*/

#include "KinKal/KTLine.hh"
#include "KinKal/BField.hh"
#include "KinKal/POCAUtil.hh"
#include "Math/AxisAngle.h"
#include "Math/VectorUtil.h"
#include <math.h>
#include <stdexcept>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
  typedef ROOT::Math::SVector<double,3> SVec3;
  vector<string> KTLine::paramTitles_ = {
      "Transverse DOCA to Z Axis (d_{0})", "Azimuth of POCA (#phi_{0})",
      "Z at POCA (z_{0})", "Cos #theta", "Time at POCA (t_{0})", "Momentum"};

  vector<string> KTLine::paramNames_ = {"d_{0}", "#phi_{0}", "z_{0}",
                                       "cos(#theta)", "t_{0}", "mom"};

  vector<string> KTLine::paramUnits_ = {"mm", "radians", "mm", "", "ns","MeV/c"};

  std::vector<std::string> const &KTLine::paramUnits() { return paramUnits_; }
  std::vector<std::string> const &KTLine::paramNames() { return paramNames_; }
  std::vector<std::string> const &KTLine::paramTitles() { return paramTitles_; }

  std::string const &KTLine::paramName(ParamIndex index) {
    return paramNames_[static_cast<size_t>(index)];
  }
  std::string const &KTLine::paramTitle(ParamIndex index) {
    return paramTitles_[static_cast<size_t>(index)];
  }
  std::string const &KTLine::paramUnit(ParamIndex index) {
    return paramUnits_[static_cast<size_t>(index)];
  }

  string KTLine::trajName_("KTLine");
  string const &KTLine::trajName() { return trajName_; }

 KTLine::KTLine( Vec4 const& pos0, Mom4 const& mom0, int charge, double bnom, TRange const& range) : KTLine(pos0,mom0,charge,Vec3(0.0,0.0,bnom),range) {}

  KTLine::KTLine(Vec4 const &pos0, Mom4 const &mom0, int charge, Vec3 const &bnom,
  TRange const &trange) :  bnom_(bnom), mass_(mom0.M()), charge_(charge), trange_(trange)
  {
    double mommag = mom0.R();
    double speed = ( mommag/ mom0.E()) * CLHEP::c_light;
    Vec3 dir = mom0.Vect().Unit();
    
    static const Vec3 zdir(0.0, 0.0, 1.0);
    static const Vec3 zpos(0.0, 0.0, 0.0);

// calculate POCA to the Z axis.  This is the reference for the parameters
    POCAUtil poca(pos0.Vect(), dir, zpos, zdir);
    double flen = dir.Dot(pos0.Vect()-poca.point1()); // flight length from reference to POCA
    Vec3 pca = poca.point1()-poca.point2(); // vector from Z axis to POCA
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
  KTLine can take in Momentum externally as a 4-vector or calculate it based. You
  can initialize the line with an origin (pos0) or the trajectory parameters
  (pdata)
  */

  KTLine::KTLine(PDATA const &pdata, KTLine const &other) : KTLine(other) {
    pars_ = pdata;
  }

  KTLine::KTLine(StateVector const& pstate, double time, double mass, int charge, Vec3 const& bnom, TRange const& range) :
    KTLine(Vec4(pstate.position().X(),pstate.position().Y(),pstate.position().Z(),time),
	Mom4(pstate.momentum().X(),pstate.momentum().Y(),pstate.momentum().Z(),mass),
	charge,bnom,range) 
  {}

  KTLine::KTLine(StateVectorMeasurement const& pstate, double time, double mass, int charge, Vec3 const& bnom, TRange const& range) :
  KTLine(pstate.stateVector(),time,mass,charge,bnom,range) {
  // derive the parameter space covariance from the global state space covariance
    DPDS dpds = dPardState(time);
    pars_.covariance() = ROOT::Math::Similarity(dpds,pstate.stateCovariance());
  }

  KTLine::KTLine(KTLine const& other, Vec3 const& bnom, double trot) : KTLine(other) {
  }

  void KTLine::position(Vec4 &pos) const {
    Vec3 pos3 = position(pos.T());
    pos.SetXYZT(pos3.X(), pos3.Y(), pos3.Z(), pos.T());
  }

  Vec3 KTLine::position(double time) const {
    return (pos0() + flightLength(time) * dir());
  }

  Vec4 KTLine::pos4(double time) const {
    Vec3 temp = position(time);
    return Vec4(temp.X(), temp.Y(), temp.Z(), time);
  }

  void KTLine::momentum(double tval, Mom4 &mom4) const {
    Vec3 dir = direction(tval);
    double momval = mom();
    mom4.SetPx(momval * dir.x());
    mom4.SetPy(momval * dir.y());
    mom4.SetPz(momval * dir.z());
    mom4.SetM(mass_);
  }

  Mom4 KTLine::momentum(double tval) const {
    Mom4 momvec;
    momentum(tval,momvec);
    return momvec;
  }
  StateVector KTLine::state(double time) const {
    return StateVector(position(time),momentum(time).Vect());
  }

  StateVectorMeasurement KTLine::measurementState(double time) const {
    // express the parameter space covariance in global state space
    DSDP dsdp = dStatedPar(time);
    return StateVectorMeasurement(state(time),ROOT::Math::Similarity(dsdp,pars_.covariance()));
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
  Vec3 KTLine::direction(double t, LocalBasis::LocDir mdir) const {

    switch (mdir) {
    case LocalBasis::perpdir: // purely polar change theta 1 = theta
      return Vec3(cosTheta() * cosPhi0(), cosTheta() * sinPhi0(), -1 * sinTheta());
      break;
    case LocalBasis::phidir: // purely transverse theta2 = -phi()*sin(theta)
      return Vec3(-sinPhi0(), cosPhi0(), 0.0);
      break;
    case LocalBasis::momdir: // along momentum
      return dir();
      break;
    default:
      throw std::invalid_argument("Invalid direction");
    }
  }

  DSDP KTLine::dPardState(double time) const{
  // aggregate state from separate X and M derivatives; parameter space is row
    KTLine::KTLine::DPDV dPdX = dPardX(time);
    KTLine::KTLine::DPDV dPdM = dPardM(time);
    DPDS dpds;
    dpds.Place_at(dPdX,0,0);
    dpds.Place_at(dPdM,0,3);
    return dpds;
  }

  DPDS KTLine::dStatedPar(double time) const {
  // aggregate state from separate X and M derivatives; parameter space is column
    KTLine::KTLine::DVDP dXdP = dXdPar(time);
    KTLine::KTLine::DVDP dMdP = dMdPar(time);
    DSDP dsdp;
    dsdp.Place_at(dXdP,0,0);
    dsdp.Place_at(dMdP,3,0);
    return dsdp;
  }

  KTLine::DVDP KTLine::dXdPar(double time) const { 
    double deltat = time-t0();
    double sinT = sinTheta();
    double cotT = 1.0/tanTheta();
    double cosT = cosTheta();
    double sinF = sin(phi0());
    double cosF = cos(phi0());
    double spd = speed();
    double gam = gamma();
    SVec3 dX_dd0(-sinF, cosF, 0.0);
    SVec3 dX_dphi0 = (sinT*spd*deltat)*dX_dd0 - d0()*SVec3(cosF,sinF,0.0);
    SVec3 dX_dz0 (0.0,0.0,1.0);
    SVec3 dX_dcost = (spd*deltat)*SVec3(-cotT*cosF,-cotT*sinF,1.0);
    SVec3 dX_dt0 = -spd*SVec3(sinT*cosF,sinT*sinF,cosT);
    SVec3 dX_dmom = -(deltat/(gam*gam*mom()))*dX_dt0;

    KTLine::KTLine::DVDP dXdP;
    dXdP.Place_in_col(dX_dd0,0,d0_);
    dXdP.Place_in_col(dX_dphi0,0,phi0_);
    dXdP.Place_in_col(dX_dz0,0,z0_);
    dXdP.Place_in_col(dX_dcost,0,cost_);
    dXdP.Place_in_col(dX_dt0,0,t0_);
    dXdP.Place_in_col(dX_dmom,0,mom_);
    return dXdP;
  }
  KTLine::DVDP KTLine::dMdPar(double time) const {
    double sinT = sinTheta();
    double cotT = 1.0/tanTheta();
    double cosT = cosTheta();
    double sinF = sin(phi0());
    double cosF = cos(phi0());

// note: dMdd0 = dMdz0 = dMdt0 = 0
    SVec3 dM_dphi0 = (mom()*sinT)*SVec3(-sinF,cosF,0);
    SVec3 dM_dcost = mom()*SVec3(-cotT*cosF,-cotT*sinF,1.0);
    SVec3 dM_dmom = SVec3(sinT*cosF,sinT*sinF,cosT);

    KTLine::KTLine::DVDP dMdP;
    dMdP.Place_in_col(dM_dphi0,0,phi0_);
    dMdP.Place_in_col(dM_dcost,0,cost_);
    dMdP.Place_in_col(dM_dmom,0,mom_);
    return dMdP;
  }

  KTLine::DPDV KTLine::dPardX(double time) const {
    double sinT = sinTheta();
    double cotT = 1.0/tanTheta();
    double sinF = sin(phi0());
    double cosF = cos(phi0());
    double E = energy();
// note: dCosTdX = dPhi0dX = dmom_dX = 0
    SVec3 dd0_dX = SVec3(-sinF,cosF,0.0);
    SVec3 dz0_dX = SVec3(-cotT*cosF,-cotT*sinF,1.0);
    SVec3 dt0_dX = -E*SVec3(cosF,sinF,0.0)/(mom()*sinT*CLHEP::c_light);
    KTLine::DPDV dPdX;
    dPdX.Place_in_row(dd0_dX,d0_,0);
    dPdX.Place_in_row(dz0_dX,z0_,0);
    dPdX.Place_in_row(dt0_dX,t0_,0);

    return dPdX;
  }

  KTLine::DPDV KTLine::dPardM(double time) const { 
    double sinT = sinTheta();
    double cosT = cosTheta();
    double cotT = 1.0/tanTheta();
    double sinF = sin(phi0());
    double cosF = cos(phi0());
    double cos2F = cosF*cosF-sinF*sinF;
    double sin2F = 2*cosF*sinF;
    Vec3 momv = momentum(time).Vect();
    Vec3 pos = position(time);
    static Vec3 zdir(0.0,0.0,1.0);
    Vec3 momt = VectorUtil::PerpVector(momv,zdir);
    double momt2 = momt.Mag2();
    double xmt = pos.Dot(momt);
    double E = energy();
    SVec3 dmom_dM(sinT*cosF, sinT*sinF, cosT);
    SVec3 dcost_dM = (1.0/mom())*(SVec3(0.0,0.0,1.0) - cosT*dmom_dM);
    SVec3 dphi0_dM = (1.0/(mom()*sinT))*SVec3(-sinF,cosF,0.0);
    SVec3 dt0_dM = (1.0/(momt2*CLHEP::c_light))*(
	xmt*( (2.0*E/momt2)*SVec3(momv.X(),momv.Y(),0.0) 
	  - (1.0/E)*( SVec3(momv.X(),momv.Y(),momv.Z())) )
	-  E*SVec3(pos.X(),pos.Y(),0.0));
    SVec3 dz0_dM = (1.0/(mom()*sinT))*(SVec3(cotT*(cos2F*pos.X() + sin2F*pos.Y()),cotT*(sin2F*pos.X()-cos2F*pos.Y()),-cosF*pos.X()-sinF*pos.Y()));
    SVec3 dd0_dM  = ( xmt/momt2 )* SVec3(sinF, -cosF, 0.0);
    KTLine::DPDV dPdM;
    dPdM.Place_in_row(dmom_dM,mom_,0);
    dPdM.Place_in_row(dcost_dM,cost_,0);
    dPdM.Place_in_row(dphi0_dM,phi0_,0);
    dPdM.Place_in_row(dz0_dM,z0_,0);
    dPdM.Place_in_row(dt0_dM,t0_,0);
    dPdM.Place_in_row(dd0_dM,d0_,0);
    return dPdM;
  }

  // derivatives of momentum projected along the given basis WRT the parameters
  KTLine::DVEC KTLine::momDeriv(double time, LocalBasis::LocDir mdir) const {
    KTLine::DPDV dPdM = dPardM(time);
    auto dir = direction(time,mdir);
    double mommag = momentumMag(time);
    return mommag*(dPdM*SVec3(dir.X(), dir.Y(), dir.Z()));
  }

  void KTLine::print(ostream &ost, int detail) const {
    auto perr = params().diagonal();
    ost << " KTLine " << range() << " parameters: ";
    for (size_t ipar = 0; ipar < KTLine::npars_; ipar++) {
      ost << KTLine::paramName(static_cast<ParamIndex>(ipar)) << " "
          << paramVal(ipar) << " +- " << perr(ipar);
      if (ipar < KTLine::npars_ - 1)
        ost << " ";
    }
    ost << " with rotation around Bnom " << bnom_ << endl;
  }

  ostream &operator<<(ostream &ost, KTLine const &lhel) {
    lhel.print(ost, 0);
    return ost;
  }

} // namespace KinKal
