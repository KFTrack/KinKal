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
    // given a drift DOCA and direction, compute mean drift time, local speed, and expected RMS of drift time, and the local drift speed
      virtual void distanceToTime(Pol2 const& drift, float& tdrift, float& tdriftrms, float& dspeed) const = 0;
      virtual float averageDriftSpeed() const = 0; // average drift speed
  };

  // simple implementation of the above using a constant drift velocity and no ExB effects.  Used for testing
  class CVD2T  : public D2T {
    public:
      virtual void distanceToTime(Pol2 const& drift, float& tdrift, float& tdriftrms, float& dspeed) const override {
	tdrift  = drift.R()/s_;
	tdriftrms = dt_; 
	dspeed = s_;
      }
      // provide seed (mm/ns) and time RMS (ns) on construction
      CVD2T(float s, float dt) :s_(s), dt_(dt) {}
      virtual float averageDriftSpeed() const override { return s_; }
    private:
      float s_; // drift speed
      float dt_; // constant drift time error
  };

}
#endif
