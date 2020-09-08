#ifndef KinKal_THit_hh
#define KinKal_THit_hh
//
//  Base class to describe a time hit (combination of time and position measurements)
//  Its interface defines how the measurement is translated into a residual (1=dimensional comparison with a trajectory)
//  The template argument is the kinematic trajectory base of the KKTrk
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/Residual.hh"
#include "KinKal/DXing.hh"
#include "KinKal/PKTraj.hh"
#include "KinKal/KKConfig.hh"
#include <memory>
#include <ostream>

namespace KinKal {
  template <class KTRAJ> class THit {
    public:
      using PKTRAJ = PKTraj<KTRAJ>;
      using DXING = DXing<KTRAJ>;
      using DXINGPTR = std::shared_ptr<DXING>;
     // default
      THit(bool active=true) : active_(active) {}
      // optionally create with an associated detector material crossing
      THit(DXINGPTR const& dxing,bool active=true) : dxing_(dxing), active_(active) {}
      virtual ~THit(){}
      // compute residual and errors WRT a predicted trajectory
      virtual void resid(PKTRAJ const& pktraj, Residual& resid) const =0;
      // count number of degrees of freedom constrained by this measurement (typically 1)
      virtual unsigned nDOF() const = 0;
      // update, and compute residual
      virtual void update(PKTRAJ const& pktraj, MIConfig const& config, Residual& resid) = 0;
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

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, THit<KTRAJ> const& thit) {
    thit.print(ost,0);
    return ost;
 }

}
#endif

