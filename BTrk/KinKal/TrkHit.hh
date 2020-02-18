#ifndef KinKal_TrkHit_hh
#define KinKal_TrkHit_hh
//
//  Base class to describe a time hit (simultaneous time and position measurement)
//  It is templated on the # of dimensions constrained, typically 1 (2 for pixels)
//  Its interface defines how the measurement is translated into a residual
//  Used as part of the kinematic Kalman fit
//
#include "BTrk/KinKal/TLine.hh"
#include "BTrk/KinKal/Context.hh"
#include "BTrk/KinKal/TPOCABase.hh"
#include "Math/SVector.h"
#include "Math/SMatrix.h"
namespace KinKal {

  template <unsigned HDIM=1> class Residual {
    public:
      typedef ROOT::Math::SVector<double,HDIM> RVec; // residual vector type
      typedef ROOT::Math::SMatrix<double,HDIM,HDIM,ROOT::Math::MatRepSym<double,HDIM> > RCov;  // associated covariance type
      // accessors
      RVec const& residual() const { return rvec_; }
      RCov const& covariance() const  { return rcov_; }
      // construct from data
      Residual(RVec const& rvec, RCov const& rcov) : rvec_(rvec), rcov_(rcov) {}
    private:
      RVec rvec_; // residual value
      RCov rcov_; // covariance matrix on residual
  };

  template <unsigned HDIM=1> class TrkHit {
    public:
      typedef ROOT::Math::SMatrix<double,HDIM,2> RDer; // type for derivative of residual WRT DOCA (0) and time (2)
      // construct from a local trajectory and the overall context
      TrkHit(TLine const& straj, Context const& context, bool active=true) : straj_(straj), context_(context), active_(active) {} 
      // The TTraj describing this measurement
      TLine const& sensorTraj() const { return straj_; }
      // Translate TPOCA into a residual and derivatives
      // Return value indicates if the position is consistent with the measurement, within the given number of sigmas
      virtual bool resid(TPOCABase const& tpoca, Residual<HDIM>& resid, RDer const& dRdDT, double nsigma) const =0;
      // Update the state given the most recent TPOCA
      virtual void update(TPOCABase const& tpoca) const = 0;
      bool active() const { return active_; }
    private:
      TrkHit() = delete;
      TLine straj_; // local linear approximation to this sensor geometry
      Context const& context_; // context,including BField for ExB effects
      bool active_; // is this hit active or not
  };
}
#endif

