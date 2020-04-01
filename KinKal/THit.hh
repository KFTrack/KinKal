#ifndef KinKal_THit_hh
#define KinKal_THit_hh
//
//  Base class to describe a time hit (simultaneous time and position measurement)
//  Its interface defines how the measurement is translated into a residual (1=dimensional comparison with a trajectory)
//  The template argument is the kinematic trajectory base of the KKTrk
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/BField.hh"
#include "KinKal/Residual.hh"
#include "KinKal/DXing.hh"
#include "KinKal/PKTraj.hh"
namespace KinKal {
  template <class KTRAJ> class THit {
    public:
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef DXing<KTRAJ> DXING;
      typedef typename KTRAJ::PDER PDER; // forward derivative type from the particle trajectory
      // construct from a local trajectory and the overall bfield
      THit(BField const& bfield, bool active=true) : bfield_(bfield), active_(active) {} 
      virtual ~THit(){}
      // compute a residual and derivatives from this measurement given the reference trajectory 
      virtual void resid(Residual& resid) const =0;
      // Update the internal state of this hit the most recent prediction
      virtual void update(PKTRAJ const& pktraj) = 0;
      // count number of degrees of freedom constrained by this measurement (typically 1)
      virtual unsigned nDOF() const = 0;
      // check the consistency of ancillary information of this measurement not used in the residual computation with the prediction
      // return value is the dimensionless number of sigma outside range, 0.0 = perfecty consistent, 1.0 is '1 sigma' tension
      virtual float tension() const = 0;
      // access to the TPOCA calculation result and derivatives
      virtual TPocaBase const& poca() const = 0;
      virtual PDER const& dDdP() const = 0;
      virtual PDER const& dTdP() const = 0;
      BField const& bfield() const { return bfield_; }
      bool isActive() const { return active_; }
      // associated material information
      virtual DXING* detCrossing() = 0; 
    private:
      THit() = delete;
      BField const& bfield_; // BField for ExB effects
      bool active_; // is this hit active or not
  };
}
#endif

