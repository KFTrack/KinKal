#ifndef KinKal_THit_hh
#define KinKal_THit_hh
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
  class TLine;
  class THit {
    public:
      enum DerivType {dDOCA=0,dt=1};
      typedef ROOT::Math::SVector<double,2> RDer; // type for derivative of residual WRT DOCA (0) and time (1)
      // construct from a local trajectory and the overall context
      THit(Context const& context, bool active=true) : context_(context), active_(active) {} 
      virtual ~THit(){}
      // The TLine describing this measurement
      virtual TLine const& sensorTraj() const = 0;
      // Translate TPOCA into a residual and derivatives
      // Return value indicates if the position is consistent with the measurement, within the given number of sigmas
      virtual bool resid(TPOCABase const& tpoca, Residual& resid, RDer const& dRdDT, double nsigma) const =0;
      // Update the state given the most recent TPOCA
      virtual void update(TPOCABase const& tpoca) const = 0;
      virtual unsigned nDOF() const = 0;
      bool active() const { return active_; }
      Context const& context() const { return context_; }
    private:
      THit() = delete;
      Context const& context_; // context,including BField for ExB effects
      bool active_; // is this hit active or not
  };
}
#endif

