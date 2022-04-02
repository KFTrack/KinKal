#ifndef KinKal_SimpleWireHit_hh
#define KinKal_SimpleWireHit_hh
//
// Simple implementation of a wire hit, for testing purpopses
//
#include "KinKal/Detector/WireHit.hh"
namespace KinKal {

  template <class KTRAJ> class SimpleWireHit : public WireHit<KTRAJ> {
    public:
      using WIREHIT = WireHit<KTRAJ>;
      using Dimension = typename WireHit<KTRAJ>::Dimension;
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using PTCA = PiecewiseClosestApproach<KTRAJ,Line>;

      SimpleWireHit(BFieldMap const& bfield, PTCA const& ptca, WireHitState const& whstate, double mindoca,
          double driftspeed, double tvar, double rcell);
      // override updating.  I have to override both since they have the same name
      void update(PKTRAJ const& pktraj) override;
      void update(PKTRAJ const& pktraj, MetaIterConfig const& config) override;
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
      // allow the updater access
      friend class SimpleWireHitUpdater;
  };

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

  // simple updater based on DOCA
  class SimpleWireHitUpdater {
    public:
      SimpleWireHitUpdater(double mindoca,double maxdoca ) : mindoca_(mindoca), maxdoca_(maxdoca) {}
      template <class KTRAJ> void update(SimpleWireHit<KTRAJ>& swh) const;
    private:
      double mindoca_; // minimum DOCA value to sign LR ambiguity
      double maxdoca_; // maximum DOCA to still use a hit
  };

  template <class KTRAJ> void SimpleWireHitUpdater::update(SimpleWireHit<KTRAJ>& swh) const {
    swh.mindoca_ = std::min(mindoca_,swh.cellRadius());
    WireHitState::State state;
    if(swh.closestApproach().usable()){
      double doca = swh.closestApproach().doca();
//      auto chisq = swh.cminprob_ hisquared(); // unbiased chisquared
//      if(fabs(doca) > maxdoca_ || chisq.probability() < minprob_ ) {
      if(fabs(doca) > maxdoca_ ) {
        state = WireHitState::inactive; // disable the hit if it's an outlier
      } else if(fabs(doca) > mindoca_ ) {
        state = doca > 0.0 ? WireHitState::right : WireHitState::left;
      } else {
        // too close to the wire: don't try to disambiguate LR sign
        state = WireHitState::null;
      }
    } else {
      state = WireHitState::inactive;
    }
    swh.setState(state);
  };

  template <class KTRAJ> SimpleWireHit<KTRAJ>::SimpleWireHit(BFieldMap const& bfield, PTCA const& ptca, WireHitState const& whstate,
      double mindoca, double driftspeed, double tvar, double rcell) :
    WIREHIT(bfield,ptca,whstate), mindoca_(mindoca), dvel_(driftspeed), tvar_(tvar), rcell_(rcell) {}

  template <class KTRAJ> void SimpleWireHit<KTRAJ>::distanceToTime(POL2 const& drift, DriftInfo& dinfo) const {
    // simply translate distance to time using the fixed velocity
    dinfo.tdrift_ = drift.R()/dvel_;
    dinfo.vdrift_ = dvel_;
    dinfo.tdriftvar_ = tvar_;
  }

  template <class KTRAJ> void SimpleWireHit<KTRAJ>::update(PKTRAJ const& pktraj) {
    WIREHIT::update(pktraj);
  }

  template <class KTRAJ> void SimpleWireHit<KTRAJ>::update(PKTRAJ const& pktraj,MetaIterConfig const& miconfig) {
  // look for an updater; if it's there, update the state
    auto swhu = miconfig.findUpdater<SimpleWireHitUpdater>();
    if(swhu != 0){
//      auto tpoca = WIREHIT::updatePTCA(pktraj);
//      WIREHIT::updateDrift(tpoca);
//      WIREHIT::update(pktraj,miconfig);
      swhu->update(*this);
    }
    WIREHIT::update(pktraj,miconfig);
  }

}
#endif
