#ifndef KinKal_SimpleWireHit_hh
#define KinKal_SimpleWireHit_hh
//
// Simple implementation of a wire hit, for testing purpopses
//
#include "KinKal/Detector/ResidualHit.hh"
#include "KinKal/Examples/DOCAWireHitUpdater.hh"
#include "KinKal/Examples/WireHitStructs.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/SensorLine.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/Trajectory/ClosestApproach.hh"
#include "KinKal/General/BFieldMap.hh"
#include <array>
#include <stdexcept>
namespace KinKal {

  template <class KTRAJ> class SimpleWireHit : public ResidualHit<KTRAJ> {
    public:
      using HIT = Hit<KTRAJ>;
      using PCA = PiecewiseClosestApproach<KTRAJ,SensorLine>;
      using CA = ClosestApproach<KTRAJ,SensorLine>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      using PTRAJ = ParticleTrajectory<KTRAJ>;
      enum Dimension { dresid=0, tresid=1, lresid=2};  // residual dimensions

      SimpleWireHit(BFieldMap const& bfield, PCA const& pca, WireHitState const& whstate, double mindoca,
          double driftspeed, double tvar, double tot, double totvar, double longval, double longvar, double rcell,int id);
      unsigned nResid() const override { return 3; } // 2 residuals
      double time() const override { return ca_.particleToca(); }
      Residual const& refResidual(unsigned ires=dresid) const override;
      void updateReference(PTRAJ const& ptraj) override;
      KTRAJPTR const& refTrajPtr() const override { return ca_.particleTrajPtr(); }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // Use dedicated updater
      void updateState(MetaIterConfig const& config,bool first) override;
      // specific to SimpleWireHit: this has a constant drift speed
      double cellRadius() const { return rcell_; }
      virtual ~SimpleWireHit(){}
      double driftVelocity() const { return dvel_; }
      double timeVariance() const { return tvar_; }
      double minDOCA() const { return mindoca_; }
      int id() const { return id_; }
      CA unbiasedClosestApproach() const;
      auto const& closestApproach() const { return ca_; }
      auto const& hitState() const { return whstate_; }
      auto const& wire() const { return wire_; }
      auto const& bfield() const { return bfield_; }
      auto precision() const { return ca_.precision(); }
    private:
      BFieldMap const& bfield_; // drift calculation requires the BField for ExB effects
      WireHitState whstate_; // current state
      SensorLine wire_; // local linear approximation to the wire of this hit, encoding all (local) position and time information.
                  // the start time is the measurement time, the direction is from
                  // the physical source of the signal (particle) to the measurement recording location (electronics), the direction magnitude
                  // is the effective signal propagation velocity along the wire, and the time range describes the active wire length
                  // (when multiplied by the propagation velocity).
      CA ca_; // reference time and position of closest approach to the wire; this is generally biased by the hit
      std::array<Residual,3> rresid_; // residuals WRT most recent reference
      double mindoca_; // effective minimum DOCA used when assigning LR ambiguity, used to define null hit properties
      double dvel_; // constant drift speed
      double tvar_; // constant time variance
      double tot_, totvar_; // TimeOverThreshold and variance
      double longval_, longvar_; // longitudinal value and variance
      double rcell_; // straw radius
      int id_; // id
      void updateResiduals();
  };

  //trivial 'updater' that sets the wire hit state to null
  class NullWireHitUpdater {
    public:
      void updateWireHitState(WireHitState &state) const { state.state_ = WireHitState::null; }
  };
  //trivial 'updater' that enables longitudinal measurement
  class LongWireHitUpdater {
    public:
      void updateWireHitState(WireHitState &state) const {state.lstate_ = WireHitState::longactive; }
  };

  template <class KTRAJ> SimpleWireHit<KTRAJ>::SimpleWireHit(BFieldMap const& bfield, PCA const& pca, WireHitState const& whstate,
      double mindoca, double driftspeed, double tvar, double tot, double totvar, double longval, double longvar, double rcell, int id) :
    bfield_(bfield),
    whstate_(whstate), wire_(pca.sensorTraj()),
    ca_(pca.localTraj(),wire_,pca.precision(),pca.tpData(),pca.dDdP(),pca.dTdP(),pca.dLdP()), // must be explicit to get the right sensor traj reference
    mindoca_(mindoca), dvel_(driftspeed), tvar_(tvar), tot_(tot), totvar_(totvar), longval_(longval), longvar_(longvar), rcell_(rcell), id_(id) {
    }

  template <class KTRAJ> void SimpleWireHit<KTRAJ>::updateReference(PTRAJ const& ptraj) {
    // if we already computed PCA in the previous iteration, use that to set the hint.  This speeds convergence
    // otherwise use the time at the center of the wire
    CAHint tphint = ca_.usable() ?  ca_.hint() : CAHint(wire_.timeAtMidpoint(),wire_.timeAtMidpoint());
    PCA pca(ptraj,wire_,tphint,precision());
    ca_ = pca.localClosestApproach();
    if(!ca_.usable())throw std::runtime_error("WireHit TPOCA failure");
  }

  template <class KTRAJ> void SimpleWireHit<KTRAJ>::updateState(MetaIterConfig const& miconfig, bool first) {
    if(first){
      // look for an updater; if found, use it to update the state
      auto nwhu = miconfig.findUpdater<NullWireHitUpdater>();
      auto twhu = miconfig.findUpdater<LongWireHitUpdater>();
      auto dwhu = miconfig.findUpdater<DOCAWireHitUpdater>();
      if(nwhu != 0 && dwhu != 0)throw std::invalid_argument(">1 SimpleWireHit updater specified");
      if(twhu != 0){
        twhu->updateWireHitState(whstate_);
      }
      if(nwhu != 0){
        mindoca_ = cellRadius();
        nwhu->updateWireHitState(whstate_);
        // set the residuals based on this state
      } else if(dwhu != 0){
        // update minDoca (for null ambiguity error estimate)
        mindoca_ = std::min(dwhu->minDOCA(),cellRadius());
        // compute the unbiased closest approach.  This is brute-force
        // a more clever solution is to linearly correct the residuals for the change in parameters
        auto uca = this->unbiasedClosestApproach();
        if (uca.usable())
          dwhu->updateWireHitState(whstate_,uca.doca());
        else
          whstate_ = WireHitState(WireHitState::inactive);
      }
    }
    rresid_[tresid] = rresid_[dresid] = rresid_[lresid] = Residual();
    if(whstate_.active()){
      rresid_[tresid] = Residual(ca_.deltaT() - tot_, totvar_,0.0,true,ca_.dTdP()); // always constrain to TOT; this stabilizes the fit
      if(whstate_.useDrift()){
        // translate PCA to residual. Use ambiguity assignment to convert drift time to a drift radius
        double dr = dvel_*whstate_.lrSign()*ca_.deltaT() -ca_.doca();
        DVEC dRdP = dvel_*whstate_.lrSign()*ca_.dTdP() -ca_.dDdP();
        rresid_[dresid] = Residual(dr,tvar_*dvel_*dvel_,0.0,true,dRdP);
      } else {
        // interpret DOCA against the wire directly as a residuals
        double nulldvar = dvel_*dvel_*(ca_.deltaT()*ca_.deltaT()+0.8);
        rresid_[dresid] = Residual(ca_.doca(),nulldvar,0.0,true,ca_.dDdP());
      }
      if(whstate_.useLong()){
        double calong = (ca_.sensorPoca().Vect() - wire_.middle()).Dot(ca_.sensorDirection());
        double lr = calong - longval_;
        rresid_[lresid] = Residual(lr,longvar_,0.0,true,ca_.dLdP());
      }
    }
    // now update the weight
    this->updateWeight(miconfig);
  }

  template <class KTRAJ> Residual const& SimpleWireHit<KTRAJ>::refResidual(unsigned ires) const {
    if(ires >lresid)throw std::invalid_argument("Invalid residual");
    return rresid_[ires];
  }

  template <class KTRAJ> ClosestApproach<KTRAJ,SensorLine> SimpleWireHit<KTRAJ>::unbiasedClosestApproach() const {
    // compute the unbiased closest approach; this is brute force, but works
    auto const& ca = this->closestApproach();
    auto uparams = HIT::unbiasedParameters();
    KTRAJ utraj(uparams,ca.particleTraj());
    return CA(utraj,this->wire(),ca.hint(),ca.precision());
  }

  template<class KTRAJ> void SimpleWireHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    ost << " WireHit state ";
    switch(whstate_.state_) {
      case WireHitState::inactive:
        ost << "inactive";
        break;
      case WireHitState::left:
        ost << "left";
        break;
      case WireHitState::right:
        ost << "right";
        break;
      case WireHitState::null: default:
        ost << "null";
        break;
    }
    switch(whstate_.lstate_) {
      case WireHitState::longinactive:
        ost << " longinactive";
        break;
      case WireHitState::longactive:
        ost << " longactive";
        break;
    }
    if(detail > 0){
      if(rresid_[tresid].active())
        ost << " Active Time Residual " << rresid_[tresid];
      if(rresid_[dresid].active())
        ost << " Active Distance Residual " << rresid_[dresid];
      if(rresid_[lresid].active())
        ost << " Active Long Residual " << rresid_[lresid];
      ost << std::endl;
    }
    if(detail > 1) {
      ost << "Approximate Propagation speed " << wire_.speed(100) << " TPOCA " << ca_.tpData() << std::endl;
    }
  }

}
#endif
