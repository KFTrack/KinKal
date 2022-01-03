#ifndef KinKal_WireHit_hh
#define KinKal_WireHit_hh
//
//  class representing a drift wire measurement.  Implemented using PTCA between the particle traj and the wire
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/Detector/ResidualHit.hh"
#include "KinKal/Detector/WireHitStructs.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/General/BFieldMap.hh"
#include <array>
#include <stdexcept>
namespace KinKal {


  template <class KTRAJ> class WireHit : public ResidualHit<KTRAJ> {
    public:
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using PTCA = PiecewiseClosestApproach<KTRAJ,Line>;
      enum Dimension { tresid=0, dresid=1};  // residual dimensions
      // Hit interface overrrides; subclass still needs to implement state change update
      unsigned nResid() const override { return 2; } // potentially 2 residuals
      bool activeRes(unsigned ires) const override;
      Residual const& residual(unsigned ires=tresid) const override;
      double time() const override { return tpdata_.particleToca(); }
      void update(PKTRAJ const& pktraj) override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // virtual interface that must be implemented by concrete WireHit subclasses
      // given a drift DOCA and direction in the cell, compute drift time and velocity
      virtual void distanceToTime(POL2 const& drift, DriftInfo& dinfo) const = 0;
      // WireHit specific functions
      ClosestApproachData const& closestApproach() const { return tpdata_; }
      WireHitState const& hitState() const { return wstate_; }
      WireHitState& hitState() { return wstate_; }
      Residual const& timeResidual() const { return rresid_[tresid]; }
      Residual const& spaceResidual() const { return rresid_[dresid]; }
      Line const& wire() const { return wire_; }
      BFieldMap const& bfield() const { return bfield_; }
      // constructor
      WireHit(BFieldMap const& bfield, PTCA const& ptca, WireHitState const&);
      virtual ~WireHit(){}
    protected:
      void setHitState(WireHitState const& newstate) { wstate_ = newstate; }
      virtual void setResiduals(PTCA const& tpoca); // compute the Residuals; TPOCA must be already calculated
      void setPrecision(double precision) { precision_ = precision; }
    private:
      BFieldMap const& bfield_; // drift calculation requires the BField for ExB effects
      Line wire_; // local linear approximation to the wire of this hit.
      // the start time is the measurement time, the direction is from
      // the physical source of the signal (particle) towards the measurement location, the vector magnitude
      // is the effective signal propagation velocity, and the range describes the active wire length
      // (when multiplied by the propagation velocity).
      WireHitState wstate_; // current state
      // caches used in processing
      ClosestApproachData tpdata_; // reference time and distance of closest approach to the wire
      std::array<Residual,2> rresid_; // residuals WRT most recent reference
      double precision_; // precision for PTCA calculation; can change during processing schedule
  };

  template <class KTRAJ> WireHit<KTRAJ>::WireHit(BFieldMap const& bfield, PTCA const& ptca, WireHitState const& wstate) :
    bfield_(bfield), wire_(ptca.sensorTraj()), wstate_(wstate), tpdata_(ptca.tpData()), precision_(ptca.precision()) {}

  template <class KTRAJ> void WireHit<KTRAJ>::update(PKTRAJ const& pktraj) {
    // compute PTCA.  Default hint is the wire middle
    CAHint tphint(wire_.range().mid(),wire_.range().mid());
    // if we already computed PTCA in the previous iteration, use that to set the hint.  This speeds convergence
    if(tpdata_.usable()) tphint = CAHint(tpdata_.particleToca(),tpdata_.sensorToca());
    // re-compute the time point of closest approache
    PTCA tpoca(pktraj,wire_,tphint,precision_);
    if(tpoca.usable()){
      tpdata_ = tpoca.tpData();
      setResiduals(tpoca);
      this->setRefParams(pktraj.nearestPiece(tpoca.particleToca()));
    } else
      throw std::runtime_error("PTCA failure");
  }

  template <class KTRAJ> bool WireHit<KTRAJ>::activeRes(unsigned ires) const {
    if(ires ==0 && wstate_.active())
      return true;
    else if(ires ==1 && wstate_.state_ == WireHitState::null)
      return true;
    else
      return false;
  }

  template <class KTRAJ> void WireHit<KTRAJ>::setResiduals(PTCA const& tpoca) {
    // compute drift parameters.  These are used even for null-ambiguity hits
    VEC3 bvec = bfield_.fieldVect(tpoca.particlePoca().Vect());
    auto pdir = bvec.Cross(wire_.direction()).Unit(); // direction perp to wire and BFieldMap
    VEC3 dvec = tpoca.delta().Vect();
    double phi = asin(double(dvec.Unit().Dot(pdir))); // azimuth around the wire WRT the BField
    POL2 drift(fabs(tpoca.doca()), phi);
    DriftInfo dinfo;
    distanceToTime(drift, dinfo);
    if(wstate_.useDrift()){
      // compute the precise drift
      // translate PTCA to residual
      // Use ambiguity to convert drift time to a time difference.   null ambiguity means ignore drift time
      double dsign = wstate_.lrSign()*tpoca.lSign(); // overall sign is the product of assigned ambiguity and doca (angular momentum) sign
      double dt = tpoca.deltaT()-dinfo.tdrift_*dsign;
      // residual is in time, so unit dependendence on time, distance dependence is the local drift velocity
      DVEC dRdP = tpoca.dDdP()*dsign/dinfo.vdrift_ - tpoca.dTdP();
      rresid_[tresid] = Residual(dt,dinfo.tdriftvar_,dRdP);
    } else {
      // interpret DOCA and TOCA against the wire directly as residuals.  We have to take the sign out of DOCA
      DVEC dRdP = -tpoca.lSign()*tpoca.dDdP();
      double dd = tpoca.doca();
      // Interpret the maximum DOCA as a variance on the distance, assuming a flat distribution
      double nulldvar = (wstate_.maxndoca_*wstate_.maxndoca_)/3.0; // range is +- doca
      rresid_[dresid] = Residual(dd,nulldvar,dRdP);
      // correct time difference for the average drift, assuming flat drift time distribution
      double dt = tpoca.deltaT() - 0.5*wstate_.maxndoca_/dinfo.vdrift_;
      // the time constraint variance is the sum of the variance from maxdoca and from the intrinsic measurement variance
      double nulltvar = (wstate_.maxndoca_*wstate_.maxndoca_)/(dinfo.vdrift_*dinfo.vdrift_*12.0) + dinfo.tdriftvar_;
      rresid_[tresid] = Residual(dt,nulltvar,-tpoca.dTdP());
      // Note there is no correlation between distance and time residuals; the former is just from the wire position, the latter from the time measurement
    }
  }

  template <class KTRAJ> Residual const& WireHit<KTRAJ>::residual(unsigned ires) const {
    if(ires >=2)throw std::invalid_argument("Invalid residual");
    return rresid_[ires];
  }

  template<class KTRAJ> void WireHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    ost << " WireHit state ";
    switch(wstate_.state_) {
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
      ost << "Propagation speed " << wire_.speed() << " TPOCA " << tpdata_ << std::endl;
    }
  }


} // KinKal namespace
#endif
