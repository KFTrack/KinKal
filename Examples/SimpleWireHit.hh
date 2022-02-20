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
      // specific to SimpleWireHit: this has a constant drift speed
      virtual ~SimpleWireHit(){}
      double timeVariance() const { return tvar_; }
    private:
      double dvel_; // constant drift speed
      double tvar_; // constant time variance
      double rcell_; // straw radius
  };

  // simple updater based on DOCA
  struct SimpleWireHitUpdater {
    SimpleWireHitUpdater(double mindoca,double maxdoca, double minprob ) : mindoca_(mindoca), maxdoca_(maxdoca), minprob_(minprob) {}
    double mindoca_; // minimum DOCA value to sign LR ambiguity
    double maxdoca_; // maximum DOCA to still use a hit
    double minprob_; // minimum residual probability to keep using a hit
  };

  template <class KTRAJ> SimpleWireHit<KTRAJ>::SimpleWireHit(BFieldMap const& bfield, PTCA const& ptca, WireHitState const& whstate,
      double mindoca, double driftspeed, double tvar, double rcell) :
    WIREHIT(bfield,ptca,whstate,mindoca), dvel_(driftspeed), tvar_(tvar), rcell_(rcell) {}

  template <class KTRAJ> void SimpleWireHit<KTRAJ>::distanceToTime(POL2 const& drift, DriftInfo& dinfo) const {
    // simply translate distance to time using the fixed velocity
    dinfo.tdrift_ = drift.R()/dvel_;
    dinfo.vdrift_ = dvel_;
    dinfo.tdriftvar_ = tvar_;
  }

  template <class KTRAJ> void SimpleWireHit<KTRAJ>::update(PKTRAJ const& pktraj) {
    return WIREHIT::update(pktraj);
  }

  template <class KTRAJ> void SimpleWireHit<KTRAJ>::update(PKTRAJ const& pktraj,MetaIterConfig const& miconfig) {
    PTCA tpoca = WIREHIT::wirePTCA(pktraj);
    if(tpoca.usable()){
      this->tpdata_ = tpoca.tpData();
      this->setRefParams(pktraj.nearestPiece(tpoca.particleToca()));
      // find the specific updater for this meta-iteration.  There should be at most 1
      auto swhu = miconfig.findUpdater<SimpleWireHitUpdater>();
      if(swhu != 0){
        // update the WireHitState
        this->mindoca_ = std::min(swhu->mindoca_,cellRadius());;
        double doca = tpoca.doca();
        auto chisq = this->chisq();
        if(fabs(doca) > swhu->maxdoca_ || chisq.probability() < swhu->minprob_ ) {
          this->wstate_.state_ = WireHitState::inactive; // disable the hit if it's an outlier
        } else if(fabs(doca) > swhu->mindoca_ ) {
          this->wstate_.state_ = doca > 0.0 ? WireHitState::right : WireHitState::left;
        } else {
          // too close to the wire: don't try to disambiguate LR sign
          this->wstate_.state_ = WireHitState::null;
        }
      }
      WIREHIT::setResiduals(tpoca);
    } else
      throw std::runtime_error("PTCA failure");
    // OK if no swhu is found, hits may be frozen this meta-iteration
  }


}
#endif
