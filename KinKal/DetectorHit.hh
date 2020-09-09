#ifndef KinKal_DetectorHit_hh
#define KinKal_DetectorHit_hh
//
//  Base class to describe a detector hit (combination of time and/or position measurements)
//  Its interface defines how the measurement is translated into a residual (1=dimensional comparison with a trajectory)
//  The template argument is the kinematic trajectory base of the KKTrk
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/Residual.hh"
#include "KinKal/DetectorXing.hh"
#include "KinKal/ParticleTrajectory.hh"
#include "KinKal/Config.hh"
#include <memory>
#include <ostream>

namespace KinKal {
  template <class KTRAJ> class DetectorHit {
    public:
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using DXING = DetectorXing<KTRAJ>;
      using DXINGPTR = std::shared_ptr<DXING>;
     // default
      DetectorHit(bool active=true) : active_(active) {}
      // optionally create with an associated detector material crossing
      DetectorHit(DXINGPTR const& dxing,bool active=true) : dxing_(dxing), active_(active) {}
      virtual ~DetectorHit(){}
      // compute residual and errors WRT a predicted trajectory
      virtual void resid(PKTRAJ const& pktraj, Residual& resid,double precision) const =0;
      // count number of degrees of freedom constrained by this measurement (typically 1)
      virtual unsigned nDOF() const = 0;
      // update, and compute residual
      virtual void update(PKTRAJ const& pktraj, MetaIterConfig const& config, Residual& resid) = 0;
      // consistency of ancillary information not used in the residual computation
      // return value is the dimensionless number of sigma outside range, 0.0 = perfectly consistent, 1.0 is '1 sigma' tension
      virtual double tension() const = 0;
      // hits may get deactivated during the fit
      bool isActive() const { return active_; }
      bool setActivity(bool newstate) { bool retval = newstate == active_; active_ = newstate; return retval; }
      // associated material information; null use_count means no material
      DXING const& detXing() const { return *dxing_; }
      DXINGPTR const& detXingPtr() const { return dxing_; }
      bool hasMaterial() const { return (bool)dxing_; }
      virtual void print(std::ostream& ost=std::cout,int detail=0) const =0;
    private:
      DXINGPTR dxing_;
      bool active_; 
  };

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, DetectorHit<KTRAJ> const& thit) {
    thit.print(ost,0);
    return ost;
 }

}
#endif

