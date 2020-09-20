#ifndef KinKal_WireCell_hh
#define KinKal_WireCell_hh
//
// class describing a wire cell.  Eventually this could have more structure
//
#include "General/Vectors.hh"
namespace KinKal {
  class WireCell {
    public:
      virtual ~WireCell() {}
      WireCell() {}
      // disallow copy and equivalence
      WireCell(WireCell const& ) = delete; 
      WireCell& operator =(WireCell const& ) = delete;
      virtual double size() const = 0; // approximate size perpendicular to the wire
      // given a drift DOCA and direction in the cell, compute most probable drift time, expected RMS of drift time, and the local drift speed
      virtual void distanceToTime(POL2 const& drift, double& tdrift, double& tdriftvar, double& dspeed) const = 0;
      virtual double averageDriftSpeed() const = 0; // average drift speed
      virtual double maximumDriftTime() const = 0; // maximum drift
      // add functions to describe material TODO
  };

  // simple implementation; a realistic implementation needs EXB effects, non-constant drift, etc. TODO
  class SimpleCell : public WireCell {
    public:
      double size() const override { return 2*rcell_; }
      void distanceToTime(POL2 const& drift, double& tdrift, double& tdriftvar, double& dspeed) const override {
      tdrift  = drift.R()/dvel_;
      tdriftvar = tvar_; 
      dspeed = dvel_;
    }
    // provide seed (mm/ns) and time RMS (ns) on construction
    SimpleCell(double driftspeed, double tvar, double rcell) :dvel_(driftspeed), tvar_(tvar), rcell_(rcell) {}
    double averageDriftSpeed() const override { return dvel_; }
    double maximumDriftTime() const override { return rcell_/dvel_; }
    virtual ~SimpleCell(){}
    double timeVariance() const { return tvar_; }
    double radius() const { return rcell_; }
    private:
    double dvel_; // constant drift speed
    double tvar_; // constant time variance
    double rcell_; // straw radius
  };

}
#endif
