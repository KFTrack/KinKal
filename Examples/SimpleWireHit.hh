#ifndef KinKal_SimpleWireHit_hh
#define KinKal_SimpleWireHit_hh
//
// Simple implementation of a wire hit, for testing purpopses
//
#include "KinKal/Detector/ResidualHit.hh"
#include "KinKal/Examples/DOCAWireHitUpdater.hh"
#include "KinKal/Examples/WireHitStructs.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/Trajectory/ClosestApproach.hh"
#include "KinKal/General/BFieldMap.hh"
#include <array>
#include <stdexcept>
namespace KinKal {

  template <class KTRAJ> class SimpleWireHit : public ResidualHit<KTRAJ> {
    public:
      using HIT = Hit<KTRAJ>;
      using PCA = PiecewiseClosestApproach<KTRAJ,Line>;
      using CA = ClosestApproach<KTRAJ,Line>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      enum Dimension { tresid=0, dresid=1};  // residual dimensions

      SimpleWireHit(BFieldMap const& bfield, PCA const& pca, WireHitState const& whstate, double mindoca,
          double driftspeed, double tvar, double rcell,int id);
      unsigned nResid() const override { return 2; } // potentially 2 residuals
      double time() const override { return ca_.particleToca(); }
      Residual const& refResidual(unsigned ires=tresid) const override;
      void updateReference(KTRAJPTR const& ktrajptr) override;
      KTRAJPTR const& refTrajPtr() const override { return ca_.particleTrajPtr(); }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // Use dedicated updater
      void updateState(MetaIterConfig const& config,bool first) override;
      double varianceScale(unsigned ires) const override;
      // specific to SimpleWireHit: this has a constant drift speed
      double cellRadius() const { return rcell_; }
      double nullVariance(Dimension dim) const;
      double nullOffset(Dimension dim) const;
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
      Line wire_; // local linear approximation to the wire of this hit, encoding all (local) position and time information.
                  // the start time is the measurement time, the direction is from
                  // the physical source of the signal (particle) to the measurement recording location (electronics), the direction magnitude
                  // is the effective signal propagation velocity along the wire, and the time range describes the active wire length
                  // (when multiplied by the propagation velocity).
      CA ca_; // reference time and position of closest approach to the wire; this is generally biased by the hit
      std::array<Residual,2> rresid_; // residuals WRT most recent reference
      double mindoca_; // effective minimum DOCA used when assigning LR ambiguity, used to define null hit properties
      double dvel_; // constant drift speed
      double tvar_; // constant time variance
      double rcell_; // straw radius
      int id_; // id
      void updateResiduals();
 };

  //trivial 'updater' that sets the wire hit state to null
  class NullWireHitUpdater {
    public:
      WireHitState wireHitState() const { return WireHitState(WireHitState::null); }
  };

  template <class KTRAJ> SimpleWireHit<KTRAJ>::SimpleWireHit(BFieldMap const& bfield, PCA const& pca, WireHitState const& whstate,
      double mindoca, double driftspeed, double tvar, double rcell, int id) :
    bfield_(bfield),
    whstate_(whstate), wire_(pca.sensorTraj()),
    ca_(pca.localTraj(),wire_,pca.precision(),pca.tpData(),pca.dDdP(),pca.dTdP()),
    mindoca_(mindoca), dvel_(driftspeed), tvar_(tvar), rcell_(rcell), id_(id) {
    }

  template <class KTRAJ> void SimpleWireHit<KTRAJ>::updateReference(KTRAJPTR const& ktrajptr) {
    // if we already computed PCA in the previous iteration, use that to set the hint.  This speeds convergence
    // otherwise use the time at the center of the wire
    CAHint tphint = ca_.usable() ?  ca_.hint() : CAHint(wire_.range().mid(),wire_.range().mid());
    ca_ = CA(ktrajptr,wire_,tphint,precision());
    if(!ca_.usable())throw std::runtime_error("WireHit TPOCA failure");
  }

  template <class KTRAJ> double SimpleWireHit<KTRAJ>::nullVariance(Dimension dim) const {
    switch (dim) {
      case dresid: default:
        return (mindoca_*mindoca_)/3.0; // doca is signed
      case tresid:
        return (mindoca_*mindoca_)/(dvel_*dvel_*12.0); // TOCA is always larger than the crossing time
    }
  }

  template <class KTRAJ> double SimpleWireHit<KTRAJ>::nullOffset(Dimension dim) const {
    switch (dim) {
      case dresid: default:
        return 0.0;
      case tresid:
        return -0.5*mindoca_/dvel_;
    }
  }

  template <class KTRAJ> void SimpleWireHit<KTRAJ>::updateState(MetaIterConfig const& miconfig, bool first) {
    if(first){
      // look for an updater; if found, use it to update the state
      auto nwhu = miconfig.findUpdater<NullWireHitUpdater>();
      auto dwhu = miconfig.findUpdater<DOCAWireHitUpdater>();
      if(nwhu != 0 && dwhu != 0)throw std::invalid_argument(">1 SimpleWireHit updater specified");
      if(nwhu != 0){
        mindoca_ = cellRadius();
        whstate_ = nwhu->wireHitState();
        // set the residuals based on this state
      } else if(dwhu != 0){
        // update minDoca (for null ambiguity error estimate)
        mindoca_ = std::min(dwhu->minDOCA(),cellRadius());
        // compute the unbiased closest approach.  This is brute-force
        // a more clever solution is to linearly correct the residuals for the change in parameters
        auto uca = this->unbiasedClosestApproach();
        whstate_ = uca.usable() ? dwhu->wireHitState(uca.doca()) : WireHitState(WireHitState::inactive);
      }
    }
    if(whstate_.active()){
     // simply translate distance to time using the fixed velocity
      double tdrift = fabs(ca_.doca())/dvel_;
      if(whstate_.useDrift()){
        // translate PCA to residual. Use ambiguity to convert drift time to a time difference.
        double dsign = whstate_.lrSign()*ca_.lSign(); // overall sign is the product of assigned ambiguity and doca (angular momentum) sign
        double dt = ca_.deltaT()-tdrift*dsign;
        // time differnce affects the residual both through the drift distance (DOCA) and the particle arrival time at the wire (TOCA)
        DVEC dRdP = ca_.dDdP()*dsign/dvel_ - ca_.dTdP();
        rresid_[tresid] = Residual(dt,tvar_,0.0,true,dRdP);
        rresid_[dresid] = Residual();
      } else {
        // interpret DOCA against the wire directly as a residuals.  We have to take the DOCA sign out of the derivatives
        DVEC dRdP = -ca_.lSign()*ca_.dDdP();
        double dd = ca_.doca() + nullOffset(dresid);
        double nulldvar = nullVariance(dresid);
        rresid_[dresid] = Residual(dd,nulldvar,0.0,true,dRdP);
        //  interpret TOCA as a residual
        double dt = ca_.deltaT() + nullOffset(tresid);
        // the time constraint variance is the sum of the variance from maxdoca and from the intrinsic measurement variance
        double nulltvar = tvar_ + nullVariance(tresid);
        rresid_[tresid] = Residual(dt,nulltvar,0.0,true,-ca_.dTdP());
        // Note there is no correlation between distance and time residuals; the former is just from the wire position, the latter from the time measurement
      }
    } else {
      rresid_[tresid] = rresid_[dresid] = Residual();
    }
 // now update the weight
    this->updateWeight(miconfig);
  }

  template <class KTRAJ> Residual const& SimpleWireHit<KTRAJ>::refResidual(unsigned ires) const {
    if(ires >dresid)throw std::invalid_argument("Invalid residual");
    return rresid_[ires];
  }

  template <class KTRAJ> double SimpleWireHit<KTRAJ>::varianceScale(unsigned ires) const {
    if(ires >dresid)throw std::invalid_argument("Invalid residual");
    if(whstate_.useDrift() && ires == tresid)
      return tvar_;
    else
      return nullVariance((Dimension)ires);
  }

 template <class KTRAJ> ClosestApproach<KTRAJ,Line> SimpleWireHit<KTRAJ>::unbiasedClosestApproach() const {
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
    if(detail > 0){
      if(rresid_[tresid].active())
        ost << " Active Time Residual " << rresid_[tresid];
      if(rresid_[dresid].active())
        ost << " Active Distance Residual " << rresid_[dresid];
      ost << std::endl;
    }
    if(detail > 1) {
      ost << "Propagation speed " << wire_.speed() << " TPOCA " << ca_.tpData() << std::endl;
    }
  }

}
#endif
