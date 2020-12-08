#ifndef KinKal_SimpleCell_hh
#define KinKal_SimpleCell_hh
//
// simple implementation of wire cell, for testing
//
#include "KinKal/Detector/WireCell.hh"
namespace KinKal {
  // simple implementation; a realistic implementation needs EXB effects, non-constant drift, etc. TODO
  class SimpleCell : public WireCell {
    public:
      double size() const override { return rcell_; }
      void distanceToTime(POL2 const& drift, double& tdrift, double& tdriftvar, double& dspeed) const override {
      tdrift  = drift.R()/dvel_;
      tdriftvar = tvar_; 
      dspeed = dvel_;
    }
    // provide seed (mm/ns) and time RMS (ns) on construction
    SimpleCell(double driftspeed, double tvar, double rcell) :dvel_(driftspeed), tvar_(tvar), rcell_(rcell) {}
    double averageDriftSpeed() const { return dvel_; }
    double maximumDriftTime() const { return rcell_/dvel_; }
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
