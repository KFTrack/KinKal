#ifndef KinKal_D2T_hh
#define KinKal_D2T_hh
//
//  Distance-to-time functor for drift in a gas
//  As a function of the distance and azimuthal angle of the POCA
//  Used as part of the kinematic Kalman fit
//
#include "BTrk/KinKal/Types.hh"
namespace KinKal {
  class D2T { 
    public:
      double distanceToTime(PV2 const& drift) const = 0;
  };

  // simple implementation of the above using a constant drift velocity.  Used for testing
  class CVD2T  : public D2T {
    public:
      double distanceToTime(PV2 const& drift) const override { return drift.R()/v_; }
      // provide velocity on construction (mm/ns)
      CVD2T(double v) :v_(v) {}
    private:
      double v_; // drift velocity magnitude
  };

}
#endif
