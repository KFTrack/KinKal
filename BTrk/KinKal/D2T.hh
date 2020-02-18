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
      virtual double distanceToTime(Pol2 const& drift) const = 0;
      virtual double speed() const = 0; // average drift speed
  };

  // simple implementation of the above using a constant drift velocity.  Used for testing
  class CVD2T  : public D2T {
    public:
      virtual double distanceToTime(Pol2 const& drift) const override { return drift.R()/s_; }
      // provide velocity on construction (mm/ns)
      CVD2T(double s) :s_(s) {}
      virtual double speed() const override { return s_; }
    private:
      double s_; // drift speed
  };

}
#endif
