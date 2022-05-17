#ifndef KinKal_WireHit_hh
#define KinKal_WireHit_hh
//
//  class representing a drift wire measurement.  Implemented using CA between the particle traj and the wire
//
#include "KinKal/Detector/ResidualHit.hh"
#include "KinKal/Detector/WireHitStructs.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/Trajectory/ClosestApproach.hh"
#include "KinKal/General/BFieldMap.hh"
#include <array>
#include <stdexcept>
namespace KinKal {


  template <class KTRAJ> class WireHit : public ResidualHit<KTRAJ> {
    public:
      using PCA = PiecewiseClosestApproach<KTRAJ,Line>;
      using CA = ClosestApproach<KTRAJ,Line>;
      using RESIDHIT = ResidualHit<KTRAJ>;
      using HIT = Hit<KTRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      enum Dimension { tresid=0, dresid=1};  // residual dimensions
      // Hit interface overrrides; subclass still needs to implement state change update
      unsigned nResid() const override { return 2; } // potentially 2 residuals
      bool activeRes(unsigned ires) const override;
      Residual const& residual(unsigned ires=tresid) const override;
      double time() const override { return tpca_.particleToca(); }
      void updateReference(KTRAJPTR const& ktrajptr) override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // virtual interface that must be implemented by concrete WireHit subclasses
      // given a drift DOCA and direction in the cell, compute drift time and velocity
      virtual void distanceToTime(POL2 const& drift, DriftInfo& dinfo) const = 0;
      // define null ambiguity hit properties
      virtual double nullVariance(Dimension dim,DriftInfo const& dinfo) const = 0;
      virtual double nullOffset(Dimension dim,DriftInfo const& dinfo) const = 0;
      // WireHit specific functions
      auto const& closestApproach() const { return tpca_; }
      auto const& hitState() const { return whstate_; }
      auto const& timeResidual() const { return rresid_[tresid]; }
      auto const& distResidual() const { return rresid_[dresid]; }
      auto const& wire() const { return wire_; }
      auto const& bfield() const { return bfield_; }
      bool hasTimeResidual() const { return whstate_ != WireHitState::inactive; }
      bool hasDistResidual() const { return whstate_ == WireHitState::null; }
      auto precision() const { return tpca_.precision(); }
      // constructor
      WireHit(BFieldMap const& bfield, PCA const& pca, WireHitState const& whs);
      virtual ~WireHit(){}
    protected:
      // allow subclasses to update the internal state or residuals
      void setWireHitState(WireHitState::State state) { whstate_.state_ = state; }
      void updateResiduals(WireHitState const& whstate);
    private:
      BFieldMap const& bfield_; // drift calculation requires the BField for ExB effects
      WireHitState whstate_; // current state
      Line wire_; // local linear approximation to the wire of this hit, encoding all (local) position and time information.
      // the start time is the measurement time, the direction is from
      // the physical source of the signal (particle) to the measurement recording location (electronics), the direction magnitude
      // is the effective signal propagation velocity along the wire, and the time range describes the active wire length
      // (when multiplied by the propagation velocity).
      CA tpca_; // reference time and position of closest approach to the wire
      std::array<Residual,2> rresid_; // residuals WRT most recent reference
  };

  template <class KTRAJ> WireHit<KTRAJ>::WireHit(BFieldMap const& bfield, PCA const& pca, WireHitState const& wstate) :
    bfield_(bfield),
    whstate_(wstate), wire_(pca.sensorTraj()),
    tpca_(pca.localTraj(),wire_,pca.precision(),pca.tpData(),pca.dDdP(),pca.dTdP()) {
      HIT::updateReference(tpca_.particleTrajPtr());
    }

  template <class KTRAJ> bool WireHit<KTRAJ>::activeRes(unsigned ires) const {
    if(ires == tresid) return hasTimeResidual();
    if(ires == dresid) return hasDistResidual();
    return false;
  }

  template <class KTRAJ> void WireHit<KTRAJ>::updateReference(KTRAJPTR const& ktrajptr) {
    // if we already computed PCA in the previous iteration, use that to set the hint.  This speeds convergence
    // otherwise use the time at the center of the wire
    CAHint tphint = tpca_.usable() ?  tpca_.hint() : CAHint(wire_.range().mid(),wire_.range().mid());
    tpca_ = CA(ktrajptr,wire_,tphint,precision());
    if(!tpca_.usable())throw std::runtime_error("WireHit TPOCA failure");
    HIT::updateReference(ktrajptr);
    // update residuals without changing state
    updateResiduals(whstate_);
  }

  template <class KTRAJ> void WireHit<KTRAJ>::updateResiduals(WireHitState const& whstate) {
    // update the state
    whstate_ = whstate;
    // compute drift parameters.  These are used even for null-ambiguity hits
    VEC3 bvec = bfield_.fieldVect(tpca_.particlePoca().Vect());
    auto pdir = bvec.Cross(wire_.direction()).Unit(); // direction perp to wire and BFieldMap
    VEC3 dvec = tpca_.delta().Vect();
    double phi = asin(double(dvec.Unit().Dot(pdir))); // azimuth around the wire WRT the BField
    POL2 drift(fabs(tpca_.doca()), phi);
    DriftInfo dinfo;
    distanceToTime(drift, dinfo);
    if(whstate_.useDrift()){
      // translate PCA to residual. Use ambiguity to convert drift time to a time difference.
      double dsign = whstate_.lrSign()*tpca_.lSign(); // overall sign is the product of assigned ambiguity and doca (angular momentum) sign
      double dt = tpca_.deltaT()-dinfo.tdrift_*dsign;
      // time differnce affects the residual both through the drift distance (DOCA) and the particle arrival time at the wire (TOCA)
      DVEC dRdP = tpca_.dDdP()*dsign/dinfo.vdrift_ - tpca_.dTdP();
      rresid_[tresid] = Residual(dt,dinfo.tdriftvar_,dRdP);
    } else {
      // interpret DOCA against the wire directly as a residuals.  We have to take the DOCA sign out of the derivatives
      DVEC dRdP = -tpca_.lSign()*tpca_.dDdP();
      double dd = tpca_.doca() + nullOffset(dresid,dinfo);
      double nulldvar = nullVariance(dresid,dinfo);
      rresid_[dresid] = Residual(dd,nulldvar,dRdP);
      //  interpret TOCA as a residual
      double dt = tpca_.deltaT() + nullOffset(tresid,dinfo);
      // the time constraint variance is the sum of the variance from maxdoca and from the intrinsic measurement variance
      double nulltvar = dinfo.tdriftvar_ + nullVariance(tresid,dinfo);
      rresid_[tresid] = Residual(dt,nulltvar,-tpca_.dTdP());
      // Note there is no correlation between distance and time residuals; the former is just from the wire position, the latter from the time measurement
    }
  }

  template <class KTRAJ> Residual const& WireHit<KTRAJ>::residual(unsigned ires) const {
    if(ires >=2)throw std::invalid_argument("Invalid residual");
    return rresid_[ires];
  }

  template<class KTRAJ> void WireHit<KTRAJ>::print(std::ostream& ost, int detail) const {
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
      if(activeRes(tresid))
        ost << " Time Residual " << rresid_[tresid];
      if(activeRes(dresid))
        ost << " Distance Residual " << rresid_[dresid];
      ost << std::endl;
    }
    if(detail > 1) {
      ost << "Propagation speed " << wire_.speed() << " TPOCA " << tpca_.tpData() << std::endl;
    }
  }


} // KinKal namespace
#endif
