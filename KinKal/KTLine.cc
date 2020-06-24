/*
  KTLine inherits from TLine but we also want it be an instance of KKTrk. We dont need all the parameters and functions to be redefined. KTLine follows from TLine.

  Original Author: S Middleton 2020

*/

#include "KinKal/KTLine.hh"
#include "KinKal/BField.hh"
#include "Math/AxisAngle.h"
#include <math.h>
#include <stdexcept>

using namespace std;
using namespace ROOT::Math;

namespace KinKal {
    /*
    KTLine can take in Momentum externally as a 4-vector or calculate it based. You can initialize the line with an origin (pos0) or the trajectory parameters (pdata)
    */ 

  KTLine::KTLine(Vec4 const& pos0, Mom4 const& mom0, int charge, double bnom, TRange const& range) :     
    KTLine(pos0,mom0,charge,Vec3(0.0,0.0,bnom),range) {  std::cout<<" KT Constructor 1 "<<endl;  
}

  KTLine::KTLine(Vec4 const& pos0, Mom4 const& mom0, int charge, Vec3 const& bnom, TRange const& range)
  : TLine(pos0.Vect(), (mom0.Vect()/mom0.E())*CLHEP::c_light, pos0.T(), range), bnom_(bnom), pos40_(pos0), mom_(mom0), charge_(charge) {
      mass_ = mom0.M();    
      std::cout<<"Constructing KT : "<<"Mom "<<this->mom_<<" Pos "<<this->pos0_<<std::endl;
      std::cout<<"Constructed with: "<<"pdata "<<pars_.parameters()<<std::endl;
  }

  KTLine::KTLine( PDATA const& pdata, double mass, int charge, double bnom, TRange const& range)
  : KTLine(pdata,mass,charge,Vec3(0.0,0.0,bnom),range){    std::cout<<" KT Constructor 3 "<<endl;  
   } 

  KTLine::KTLine( PDATA const& pdata, double mass, int charge, Vec3 const& bnom, TRange const& range)
  : KTLine(pdata.parameters(),pdata.covariance(),mass,charge,bnom,range) {std::cout<<" KT Constructor 4 "<<endl;  
   }
  
  KTLine::KTLine(PDATA::DVEC const &pvec, PDATA::DMAT const &pcov, double mass, int charge, Vec3 const &bnom, TRange const &trange) : 
    TLine(pvec, pcov),  bnom_(bnom), mass_(mass), charge_(charge) {std::cout<<" KT Constructor 5 "<<endl;  
}

KTLine::KTLine(PDATA const& pdata, KTLine const& other) : KTLine(other) {  
    pars_ = pdata;
    //TODO  in this constructor you save the pos0 and mom of the "other" yet you have aprameters from pdata - is this consistant??? Surely you will end up with unphysical results. 
    std::cout<<"KT Constructor 6 "<<"pdata "<<pdata.parameters()<<std::endl;
  }

  string KTLine::trajName_("KTLine");  
  string const& KTLine::trajName() { return trajName_; }

  Vec4 KTLine::pos4(double time) const {
    Vec3 temp = position(time);
    return Vec4(temp.X(),temp.Y(),temp.Z(),time);
  }

  void KTLine::momentum(double tval, Mom4& mom) const{
    Vec3 dir = direction(tval);
    mom.SetPx(momentumMag(tval)*dir.x());
    mom.SetPy(momentumMag(tval)*dir.y());
    mom.SetPz(momentumMag(tval)*dir.z());
    mom.SetM(mass_);

  }

  Mom4 KTLine::momentum(double tval) const{
    Mom4 mom;
    Vec3 dir = direction(tval);
    mom.SetPx(momentumMag(tval)*dir.x());
    mom.SetPy(momentumMag(tval)*dir.y());
    mom.SetPz(momentumMag(tval)*dir.z());
    mom.SetM(mass_);
    return mom_;
  }

/*

The effects for changes in 2 perpendicular directions (theta1 = theta and
theta2 = phi()*sin(theta) can sometimes be added, as scattering in these
are uncorrelated. These axes are track specific. as cosmics are not always
coming along the same track direction it is necessary to have difference
parameterization than that used for the helix case.

alt dir = a test with the "BTrk parameterization" - just changes signs due to swithc in cos<->sin
*/
  Vec3  KTLine::direction(double t, LocalBasis::LocDir mdir) const {

    switch ( mdir ) {
    case LocalBasis::perpdir: // purely polar change theta 1 = theta
      return l2g_(Vec3(cosTheta()*sinPhi0(),cosTheta()*cosPhi0(),-1*sinTheta()));
// alt dir l2g_(Vec3(cosTheta()*cosPhi0(),cosTheta()*sinPhi0(),-1*sinTheta()));
    break;
      case LocalBasis::phidir: // purely transverse theta2 = -phi()*sin(theta)
      return l2g_(Vec3(-cosPhi0(),sinPhi0(),0.0));
//alt dir = l2g_(Vec3(-sinPhi0(),cosPhi0(),0.0));//;
    break;
      case LocalBasis::momdir: // along momentum: check.
      return l2g_(Vec3(dir().x(),dir().y(),dir().z()));
    break;
      default:
      throw std::invalid_argument("Invalid direction");
    }
  }

// derivatives of momentum projected along the given basis WRT the 5 parameters
   KTLine::DVEC KTLine::momDeriv(double time, LocalBasis::LocDir mdir) const {
    // compute some useful quantities
    double vz = CLHEP::c_light*mom().z()/mom().E();
    double l = speed()*time;
    KTLine::DVEC pder;
    cout<<"Mom deriv start params "<<pder<<endl;
    // cases
    switch ( mdir ) {
      case LocalBasis::perpdir:
	      // polar bending: change in Theta
	      pder[cost_] = -sinTheta();
	      pder[d0_] = 0;
	      pder[phi0_] = 0;
	      pder[z0_] = l*cosTheta();//alt dir =-l*cosTheta();
	      pder[t0_] = pder[z0_]/vz;
         cout<<"Mom deriv perpdir params "<<pder<<endl;
	      break;
      case LocalBasis::phidir:
	      // change in dP/dtheta1 = dP/dphi0*(-1/sintheta)
	      pder[cost_] = 0;
	      pder[d0_] = l;//alt dir = -l;
	      pder[phi0_] = -1/sinTheta();//alt dir = -1/sinTheta();
	      pder[z0_] = d0()/(sinTheta()*tanTheta());//alt dir = -d0()/(sinTheta()*tanTheta());
	      pder[t0_] = pder[z0_] /vz;
         cout<<"Mom deriv phidir params "<<pder<<endl;
	      break;
      case LocalBasis::momdir:
	      pder[cost_] = 0;
	      pder[d0_] = 0;
	      pder[phi0_] = 0;
	      pder[z0_] = 0;
	      pder[t0_] = pder[z0_] /vz;
         cout<<"Mom deriv momdir params "<<pder<<endl;
	    break;
      
      default:
        throw std::invalid_argument("Invalid direction");
        }
      return pder;
  }
  

   void KTLine::print(ostream& ost, int detail) const {
    auto perr = params().diagonal(); 
    ost << " KTLine " << range() << " parameters: ";
    for(size_t ipar=0;ipar < KTLine::npars_;ipar++){
      ost << KTLine::paramName(static_cast<KTLine::ParamIndex>(ipar) ) << " " << paramVal(ipar) << " +- " << perr(ipar);
      if(ipar < KTLine::npars_-1) ost << " ";
    }
    ost << " with rotation around Bnom " << bnom_ << endl;
  }

  ostream& operator <<(ostream& ost, KTLine const& lhel) {
    lhel.print(ost,0);
    return ost;
  }

} // KinKal namespace
