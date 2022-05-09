#ifndef KinKal_SimpleWireHit_hh
#define KinKal_SimpleWireHit_hh
//
// Simple implementation of a wire hit, for testing purpopses
//
#include "KinKal/Detector/WireHit.hh"
#include "KinKal/Examples/DOCAWireHitUpdater.hh"
namespace KinKal {

  template <class KTRAJ> class SimpleWireHit : public WireHit<KTRAJ> {
    public:
      using HIT = Hit<KTRAJ>;
      using WIREHIT = WireHit<KTRAJ>;
      using Dimension = typename WireHit<KTRAJ>::Dimension;
      using PCA = PiecewiseClosestApproach<KTRAJ,Line>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;

      SimpleWireHit(BFieldMap const& bfield, PCA const& pca, WireHitState const& whstate, double mindoca,
          double driftspeed, double tvar, double rcell);
      // Use dedicated updater
      void updateState(MetaIterConfig const& config) override;
      // WireHit interface implementations
      void distanceToTime(POL2 const& drift, DriftInfo& dinfo) const override;
      double cellRadius() const { return rcell_; }
      double nullVariance(Dimension dim,DriftInfo const& dinfo) const override;
      double nullOffset(Dimension dim,DriftInfo const& dinfo) const override;
      // specific to SimpleWireHit: this has a constant drift speed
      virtual ~SimpleWireHit(){}
      double driftVelocity() const { return dvel_; }
      double timeVariance() const { return tvar_; }
      double minDOCA() const { return mindoca_; }
    private:
      double mindoca_; // effective minimum DOCA used when assigning LR ambiguity, used to define null hit properties
      double dvel_; // constant drift speed
      double tvar_; // constant time variance
      double rcell_; // straw radius
  };

  template <class KTRAJ> SimpleWireHit<KTRAJ>::SimpleWireHit(BFieldMap const& bfield, PCA const& pca, WireHitState const& whstate,
      double mindoca, double driftspeed, double tvar, double rcell) :
    WIREHIT(bfield,pca,whstate), mindoca_(mindoca), dvel_(driftspeed), tvar_(tvar), rcell_(rcell) {}

  template <class KTRAJ> void SimpleWireHit<KTRAJ>::distanceToTime(POL2 const& drift, DriftInfo& dinfo) const {
    // simply translate distance to time using the fixed velocity
    dinfo.tdrift_ = drift.R()/dvel_;
    dinfo.vdrift_ = dvel_;
    dinfo.tdriftvar_ = tvar_;
  }

  template <class KTRAJ> double SimpleWireHit<KTRAJ>::nullVariance(Dimension dim,DriftInfo const& dinfo) const {
    switch (dim) {
      case WIREHIT::dresid: default:
        return (mindoca_*mindoca_)/3.0; // doca is signed
      case WIREHIT::tresid:
        return (mindoca_*mindoca_)/(dinfo.vdrift_*dinfo.vdrift_*12.0); // TOCA is always larger than the crossing time
    }
  }

  template <class KTRAJ> double SimpleWireHit<KTRAJ>::nullOffset(Dimension dim,DriftInfo const& dinfo) const {
    switch (dim) {
      case WIREHIT::dresid: default:
        return 0.0; // not sure if there's a better answer
      case WIREHIT::tresid:
        return -0.5*mindoca_/dinfo.vdrift_;
    }
  }

  template <class KTRAJ> void SimpleWireHit<KTRAJ>::updateState(MetaIterConfig const& miconfig) {
    // look for an updater; if it's there, update the state
    auto dwhu = miconfig.findUpdater<DOCAWireHitUpdater>();
    if(dwhu != 0){
      // update minDoca (for null ambiguity error estimate)
      mindoca_ = std::min(mindoca_,cellRadius());
      // compute the unbiased CA FIXME
      WireHitState whstate;
      if(this->closestApproach().usable()) {
        whstate = dwhu->wireHitState(this->closestApproach().doca());
      } else {
        whstate = WireHitState::inactive;
      }
      // set the residuals based on this state
      this->updateResiduals(whstate);
    }
    // update the temp.
    HIT::updateState(miconfig);
  }
}
#endif
