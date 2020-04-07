#ifndef KinKal_D2T_hh
#define KinKal_D2T_hh
//
//  Distance-to-time functor for drift in a gas
//  As a function of the distance and azimuthal angle of the POCA
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/Vectors.hh"
namespace KinKal {
  class D2T { 
    public:
    // given a drift DOCA and direction, compute mean drift time, local speed, and expected RMS of drift time, and the local drift speed
      virtual void distanceToTime(Pol2 const& drift, float& tdrift, float& tdriftvar, float& dspeed) const = 0;
      virtual float averageDriftSpeed() const = 0; // average drift speed
      virtual float maximumDriftTime() const = 0; // maximum drift
      virtual ~D2T(){}
  };

  // simple implementation of the above using a constant drift velocity and no ExB effects.  Used for testing
  class CVD2T  : public D2T {
    public:
      virtual void distanceToTime(Pol2 const& drift, float& tdrift, float& tdriftvar, float& dspeed) const override {
	tdrift  = drift.R()/dvel_;
	tdriftvar = tvar_; 
	dspeed = dvel_;
      }
      // provide seed (mm/ns) and time RMS (ns) on construction
      CVD2T(float s, float tvar, float rcell) :dvel_(s), tvar_(tvar), rcell_(rcell) {}
      virtual float averageDriftSpeed() const override { return dvel_; }
      virtual float maximumDriftTime() const override { return rcell_/dvel_; }
      virtual ~CVD2T(){}
    private:
      float dvel_; // constant drift speed
      float tvar_; // constant time variance
      float rcell_; // cell radius, used to compute maximum drift
  };

}
#endif
