#ifndef KinKal_KTLine_hh
#define KinKal_KTLine_hh
/*
  Class to join the TLine with KInter for momentum. Create fit using KTrk.
  But before that we need to add in momentum interface here.
  As KTLine inherits from TLine, there is no need to define
  position and direction.
  Original Author: S Middleton 2020
*/
#include "KinKal/BField.hh"
#include "KinKal/Vectors.hh"
#include "KinKal/TRange.hh"
#include "KinKal/PData.hh"
#include "KinKal/TLine.hh"
#include "KinKal/LocalBasis.hh"
#include "CLHEP/Units/PhysicalConstants.h"
#include "Math/Rotation3D.h"
#include <vector>
#include <stdexcept>
namespace KinKal {

  class KTLine :  public TLine {
    public:

      constexpr static ParamIndex t0Index() { return t0_; }
      typedef ROOT::Math::SVector<double,npars_> PDER; // derivative of parameters type
      typedef typename PDATA::DVEC DVEC; // derivative of parameters type
      static std::string const& trajName();
      // This also requires the nominal BField, which can be a vector (3d) or a scalar (B along z)
      KTLine(Vec4 const& pos, Mom4 const& mom, int charge, Vec3 const& bnom, TRange const& range=TRange());
      KTLine(Vec4 const& pos, Mom4 const& mom, int charge, double bnom, TRange const& range=TRange());

      KTLine(PDATA const& pdata, double mass, int charge, Vec3 const& bnom, TRange const& range=TRange());
      KTLine(PDATA const& pdata, double mass, int charge, double bnom, TRange const& range=TRange());

      KTLine(PDATA::DVEC const &pvec, PDATA::DMAT const &pcov, double mass, int charge, Vec3 const &bnom, TRange const &range = TRange());
  
      
      KTLine(PDATA const& pdata, const KTLine& ktline); 

      virtual ~KTLine() {}

      // particle momentum as a function of time
      Mom4 momentum(double time) const;
      void momentum(double t, Mom4& mom) const;

      // scalar momentum and energy in MeV/c units --> Needed for KKTrk:
      double momentumMag(double time) const  { return  gamma()*mass_*beta(); }//in MeV/c
      Mom4 mom() const { return mom_;}
      double momentumVar(double time) const  { return -1.0; }
      double energy(double time){ return mom_.E();}

      // speed in mm/ns
      void print(std::ostream& ost, int detail) const;
      void rangeInTolerance(TRange& range, BField const& bfield, double tol) const {};//infinity for striaght line

      // local momentum direction basis
      Vec3 direction(double time, LocalBasis::LocDir mdir=LocalBasis::momdir) const;
      Vec4 pos4(double time) const;
      // momentum change derivatives; this is required to instantiate a KalTrk 
      KTLine::DVEC momDeriv(double time, LocalBasis::LocDir mdir) const ;

      //some possibly useful equations:
      double mass() const{ 
        return mass_;
      }
      double ztime(double zpos) const { 
        return (t0() + zpos/((speed()*dir()).z())); 
      } //time to travel Z
      int charge() const { return charge_; }
      double beta() const { 
        return (speed()/CLHEP::c_light);
      }
      double gamma() const {
        return (1/sqrt(1-((speed()/CLHEP::c_light)*(speed()/CLHEP::c_light))));
      }
      double betaGamma() const{ 
        return beta()*gamma();
      }
      double energyBG(double time) const  { return (sqrt(mass_*mass_ + betaGamma()*betaGamma()* mass_*mass_)); }//in MeV 
      Vec3 const& bnom() const { return bnom_; }

      void invertCT(){
        charge_ *=-1;
        pars_.parameters()[t0_] *= -1.0;
      }

    private :
      Vec3 bnom_; //should be 0,0,0
      bool needsrot_; // logical flag if Bnom is parallel to global Z or not
      ROOT::Math::Rotation3D brot_; // rotation from the internal coordinate system (along B) to the global
      Mom4 pos40_, mom_; // 4 momentum vector - px,py,pz,m
      double mass_; //mass in MeV/c2
      int charge_;
      PDATA pars_;
      static std::string trajName_;
 };

}
#endif
