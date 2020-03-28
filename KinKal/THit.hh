#ifndef KinKal_THit_hh
#define KinKal_THit_hh
//
//  Base class to describe a time hit (simultaneous time and position measurement)
//  Its interface defines how the measurement is translated into a residual (1=dimensional comparison with a trajectory)
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/BField.hh"
#include "KinKal/TPocaBase.hh"
#include "KinKal/Residual.hh"
#include "KinKal/DMat.hh"
namespace KinKal {
  class TLine;
  class THit {
    public:
      // construct from a local trajectory and the overall bfield
      THit(BField const& bfield, bool active=true) : bfield_(bfield), active_(active) {} 
      virtual ~THit(){}
      // The trajectory describing this measurement's sensor.
      virtual TLine const& sensorTraj() const = 0;
      // Translate TPoca into a residual
      virtual void resid(TPocaBase const& tpoca, Residual& resid) const =0;
      // Update the state given the most recent TPoca
      virtual void update(TPocaBase const& tpoca) const = 0;
      // count number of degrees of freedom constrained by this measurement (typically 1)
      virtual unsigned nDOF() const = 0;
      // check if this POCA is within the sensor range, geometrically and temporally.
      // return value is the dimensionless number of sigma outside range, 0.0 = consistent with the range
      virtual float inRange(TPocaBase const& tpoca) const = 0;
      bool isActive() const { return active_; }
      BField const& bfield() const { return bfield_; }
      virtual const DMat* material() const { return 0; } // return the material associated with this hit
    private:
      THit() = delete;
      BField const& bfield_; // BField for ExB effects
      bool active_; // is this hit active or not
  };
}
#endif

