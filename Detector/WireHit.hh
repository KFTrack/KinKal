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
      void update(PKTRAJ const& pktraj, MetaIterConfig const& config) override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // virtual interface that must be implemented by concrete WireHit subclasses
      // given a drift DOCA and direction in the cell, compute drift time and velocity
      virtual void distanceToTime(POL2 const& drift, DriftInfo& dinfo) const = 0;
      // define null ambiguity hit properties
      virtual double nullVariance(Dimension dim,DriftInfo const& dinfo) const = 0;
      virtual double nullOffset(Dimension dim,DriftInfo const& dinfo) const = 0;
      // WireHit specific functions
      ClosestApproachData const& closestApproach() const { return tpdata_; }
      WireHitState const& hitState() const { return whstate_; }
      Residual const& timeResidual() const { return rresid_[tresid]; }
      Residual const& spaceResidual() const { return rresid_[dresid]; }
      Line const& wire() const { return wire_; }
      BFieldMap const& bfield() const { return bfield_; }
      double precision() const { return precision_; }
      // constructor
      WireHit(BFieldMap const& bfield, PTCA const& ptca, WireHitState const& whs);
      virtual ~WireHit(){}
    protected:
      void setState(WireHitState::State state) { whstate_.state_ = state; }
    private:
      WireHitState whstate_; // current state
      ClosestApproachData tpdata_; // reference time and distance of closest approach to the wire
      BFieldMap const& bfield_; // drift calculation requires the BField for ExB effects
      Line wire_; // local linear approximation to the wire of this hit.
      // the start time is the measurement time, the direction is from
      // the physical source of the signal (particle) towards the measurement location, the vector magnitude
      // is the effective signal propagation velocity, and the range describes the active wire length
      // (when multiplied by the propagation velocity).
      std::array<Residual,2> rresid_; // residuals WRT most recent reference
      double precision_; // precision for PTCA calculation; can change during processing schedule
  };

  template <class KTRAJ> WireHit<KTRAJ>::WireHit(BFieldMap const& bfield, PTCA const& ptca, WireHitState const& wstate) :
    whstate_(wstate), tpdata_(ptca.tpData()), bfield_(bfield), wire_(ptca.sensorTraj()), precision_(ptca.precision()) {}

  template <class KTRAJ> void WireHit<KTRAJ>::update(PKTRAJ const& pktraj,MetaIterConfig const& miconfig) {
    update(pktraj);
  }

  template <class KTRAJ> bool WireHit<KTRAJ>::activeRes(unsigned ires) const {
    if(ires ==0 && whstate_.active())
      return true;
    else if(ires ==1 && whstate_.state_ == WireHitState::null)
      return true;
    else
      return false;
  }

  template <class KTRAJ> void WireHit<KTRAJ>::update(PKTRAJ const& pktraj) {
    CAHint tphint(wire_.range().mid(),wire_.range().mid());
    // if we already computed PTCA in the previous iteration, use that to set the hint.  This speeds convergence
    if(tpdata_.usable()) tphint = CAHint(tpdata_.particleToca(),tpdata_.sensorToca());
    PTCA tpoca(pktraj,wire_,tphint,precision_);
    if(!tpoca.usable())throw std::runtime_error("PTCA failure");
    tpdata_ = tpoca.tpData();
    this->setRefParams(pktraj.nearestPiece(tpoca.particleToca()));
    // compute drift parameters.  These are used even for null-ambiguity hits
    VEC3 bvec = bfield_.fieldVect(tpoca.particlePoca().Vect());
    auto pdir = bvec.Cross(wire_.direction()).Unit(); // direction perp to wire and BFieldMap
    VEC3 dvec = tpoca.delta().Vect();
    double phi = asin(double(dvec.Unit().Dot(pdir))); // azimuth around the wire WRT the BField
    POL2 drift(fabs(tpoca.doca()), phi);
    DriftInfo dinfo;
    distanceToTime(drift, dinfo);
    if(whstate_.useDrift()){
      // translate PTCA to residual. Use ambiguity to convert drift time to a time difference.
      double dsign = whstate_.lrSign()*tpoca.lSign(); // overall sign is the product of assigned ambiguity and doca (angular momentum) sign
      double dt = tpoca.deltaT()-dinfo.tdrift_*dsign;
      DVEC dRdP = tpoca.dDdP()*dsign/dinfo.vdrift_ - tpoca.dTdP();
      rresid_[tresid] = Residual(dt,dinfo.tdriftvar_,dRdP);
    } else {
      // interpret DOCA against the wire directly as a residuals.  We have to take the DOCA sign out of the derivatives
      DVEC dRdP = -tpoca.lSign()*tpoca.dDdP();
      double dd = tpoca.doca() + nullOffset(dresid,dinfo);
      double nulldvar = nullVariance(dresid,dinfo);
      rresid_[dresid] = Residual(dd,nulldvar,dRdP);
      //  interpret TOCA as a residual
      double dt = tpoca.deltaT() + nullOffset(tresid,dinfo);
      // the time constraint variance is the sum of the variance from maxdoca and from the intrinsic measurement variance
      double nulltvar = dinfo.tdriftvar_ + nullVariance(tresid,dinfo);
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
      ost << "Propagation speed " << wire_.speed() << " TPOCA " << tpdata_ << std::endl;
    }
  }


} // KinKal namespace
#endif
