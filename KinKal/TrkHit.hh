#ifndef KinKal_TrkHit_hh
#define KinKal_TrkHit_hh
//
//  Base class to describe a time hit (simultaneous time and position measurement)
//  It is templated on the # of dimensions constrained, typically 1 (2 for pixels)
//  Its interface defines how the measurement is translated into a residual
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/Context.hh"
#include "KinKal/TPOCABase.hh"
#include "KinKal/Residual.hh"
namespace KinKal {
  class TTraj;
  template <unsigned HDIM=1> class TrkHit {
    public:
      typedef ROOT::Math::SMatrix<double,HDIM,2> RDer; // type for derivative of residual WRT DOCA (0) and time (1)
      // construct from a local trajectory and the overall context
      TrkHit(Context const& context, bool active=true) : context_(context), active_(active) {} 
      virtual ~TrkHit(){}
      // The TTraj describing this measurement
      virtual TTraj const& sensorTraj() const = 0;
      // Translate TPOCA into a residual and derivatives
      // Return value indicates if the position is consistent with the measurement, within the given number of sigmas
      virtual bool resid(TPOCABase const& tpoca, Residual<HDIM>& resid, RDer const& dRdDT, double nsigma) const =0;
      // Update the state given the most recent TPOCA
      virtual void update(TPOCABase const& tpoca) const = 0;
      bool active() const { return active_; }
      unsigned nDOF() const { return active() ? HDIM : 0; }
      Context const& context() const { return context_; }
    private:
      TrkHit() = delete;
      Context const& context_; // context,including BField for ExB effects
      bool active_; // is this hit active or not
  };
}
#endif

