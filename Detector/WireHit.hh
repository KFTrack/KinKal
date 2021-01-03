#ifndef KinKal_WireHit_hh
#define KinKal_WireHit_hh
//
//  class representing a drift wire measurement.  Implemented using PTCA between the particle traj and the wire
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/Detector/Hit.hh"
#include "KinKal/Detector/WireCell.hh"
#include "KinKal/Detector/ElementXing.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/General/LRAmbig.hh"
#include "KinKal/Detector/BFieldMap.hh"
#include <array>
#include <stdexcept>
namespace KinKal {

// struct describing local drift info
  struct DriftInfo {
    DriftInfo() : tdrift_(0.0), tdriftval_(0.0), dspeed_(0.0) {}
    double tdrift_; // drift time
    double tdriftvar_; // variance on drift time
    double dspeed_; // instantanious drift speed
  };

  // struct describing wire hit internal state
  struct WireHitState {
    enum LRAmbig { left=-1, null=0, right=1};  // ambiguity
    enum Dimension {none=-1, time=0, distance=1,  both=2}; // what gets constrained
    LRAmbig lrambig_; // left-right ambiguity
    Dimension dimension_; // physical dimensions being constrained
    double nullvar_; // spatial variance for null ambiguity hits
    WireHitState(LRAmbig lrambig, Dimension dim,double nvar) : lrambig_(lrambig), dimension_(dim), nullvar_(nvar) {
      if(dimension_ > time && (lrambig_ != null) throw std::invalid_argument("Inconsistant wire hit state");
    }
    WireHitState() : WireHitState(null,none,1.0) {}
  }

  template <class KTRAJ> class WireHit : public Hit<KTRAJ> {
    public:
      using HIT = Hit<KTRAJ>;
      using reftraj_ = HIT::reftraj_;
      using weight_ = HIT::weight_;
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using PTCA = PiecewiseClosestApproach<KTRAJ,Line>;
      using EXING = ElementXing<KTRAJ>;
      using EXINGPTR = std::shared_ptr<EXING>;
      
     // Hit interface overrrides; note that "update with state change" is not overridden, that must be implemented by a subclass
      unsigned nResid() const override { return 2; } // potentially 2 residuals
      bool active(unsigned ires) const override;
      Residual residual(unsigned ires=0) const override;
      double time() const override { return tpoca_.particleToca(); }
      void update(PKTRAJ const& pktraj) override;
      EXINGPTR const& detXingPtr() const override { return dxing_; }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // virtual interface that must be implemented by concrete WireHit subclasses
      // given a drift DOCA and direction in the cell, compute drift time and velocity
      virtual void distanceToTime(POL2 const& drift, DriftInfo& dinfo) const = 0;
      // WireHit specific functions
      PTCA const& closestApproach() const { return tpoca_; }
      WireHitState const& hitState() const { return wstate_; }
      Residual const& timeResidual() const { return rresid_[WireHitState::time]; }
      Residual const& spaceResidual() const { return rresid_[WireHitState::distance]; }
      Line const& wire() const { return wire_; }
      // constructor
      WireHit(BFieldMap const& bfield, Line const& wire, EXINGPTR const& dxing, WireHitState const&);
      virtual ~WireHit(){}
    protected:
      void setHitState(WireHitState const& newstate) { wstate_ = newstate; }
      virtual void setResiduals(); // compute the Residuals; TPOCA must be already calculated
    private:
      BFieldMap const& bfield_; // drift calculation requires the BField for ExB effects
      Line wire_; // local linear approximation to the wire of this hit.  The range describes the active wire length
      WireHitState wstate_; // current state
      bool active_; // active or not (pat. rec. tool)
      EXINGPTR dxing_; // material xing
      // caches used in processing
      PTCA tpoca_; // reference time and distance of closest approach to the wire
      std::array<Residual,2> rresid_; // residuals WRT most recent reference
      double precision_; // precision for PTCA calculation; can change during processing schedule
  };

  template <class KTRAJ> WireHit<KTRAJ>::WireHit(BFieldMap const& bfield, Line const& wire, EXINGPTR const& dxing, WireHitState const& wstate) : 
    bfield_(bfield), wire_(wire), dxing_(dxing), wstate_(wstate), precision_(1e-6) {}

  template <class KTRAJ> void WireHit<KTRAJ>::update(PKTRAJ const& pktraj) {
    // compute PTCA.  Default hint is the wire middle
    CAHint tphint(wire_.range().mid(),wire_.range().mid());
    // if we already computed PTCA in the previous iteration, use that to set the hint.  This speeds convergence
    if(tpoca_.usable()) tphint = CAHint(tpoca_.particleToca(),tpoca_.sensorToca());
    // re-compute the time point of closest approache
    tpoca_ = PTCA(pktraj,wire_,tphint,precision_);
    if(tpoca_.usable()){
      setResiduals();
      reftraj_ = pktraj.nearestPiece(tpoca.particleToca());
    } else
      throw std::runtime_error("PTCA failure");
  }

  template <class KTRAJ> bool WireHit<KTRAJ>::active(unsigned ires) const {
    if(ires ==0 && (wstate_.dimension_ == WireHitState::time && wstate_.dimension_ == WireHitState::both))
      return true;
    else if(ires ==1 && (wstate_.dimension == WireHitState::space)
      return true;
    else
      return false;
  }

  template <class KTRAJ> void WireHit<KTRAJ>::setResiduals() {
  // if we're using drift, convert DOCA into time
    if(wstate_.lrambig_ != WireHitState::null){
      // compute the precise drift
      // translate PTCA to residual
      VEC3 bvec = bfield_.fieldVect(tpoca_.particlePoca().Vect());
      auto pdir = bvec.Cross(wire_.dir()).Unit(); // direction perp to wire and BFieldMap
      VEC3 dvec = tpoca_.delta().Vect();
      double phi = asin(double(dvec.Unit().Dot(pdir)));
      // must use absolute DOCA to call distanceToTime
      POL2 drift(fabs(tpoca_.doca()), phi);
      DriftInfo dinfo;
      distanceToTime(drift, dinfo);
      // Use ambiguity to convert drift time to a time difference.   null ambiguity means ignore drift time
      double dsign = wstate_.lrambig_*tpoca_.lSign(); // overall sign is the product of ambiguity and doca sign
      double dt = tpoca_.deltaT()-dinfo.tdrift_*dsign;
      // residual is in time, so unit dependendence on time, distance dependence is the local drift velocity
      DVEC dRdP = tpoca_.dDdP()*dsign/dinfo.vdrift_ - tpoca_.dTdP(); 
      rresid_[WireHitState::time] = Residual(dt,dinfo.tdvar_,dRdP);
    } else {
      // interpret DOCA against the wire directly as the residual.  Sign by the angular momentum
      DVEC dRdP = tpoca_.dDdP();
      double dd = -fabs(tpoca_.doca())*tpoca_.lSign();
      rresid_[WireHitState::distance] = Residual(dd,wstate_.nullvar_,dRdP);
      if(wstate_.dimension_ == WireHitState::both){
	// add an absolute time constraint for the null ambig hits.
	// correct time difference for the average drift.  This will need to be calibrated for a real ambig resolver TODO!
	double dt = tpoca_.deltaT() - 0.5*mint;
	// convert spatial variance into temporal using the drift velocity at the wire
	POL2 drift(0.0,0.0);
	DriftInfo dinfo;
	distanceToTime(drift, dinfo);
	double dtvar = wstate_.nullvar_/(dinfo.vdrift_*dinfo.vdrift_); 
	rresid_[WireHitState::time] = Residual(dt,dtvar,-tpoca_.dTdP());
      }
    }
  }

  template <class KTRAJ> Residual const& WireHit<KTRAJ>::residual(unsigned ires) const {
    if(ires >=2)throw std::invalid_argument("Invalid residual");
    return rresid_[ires];
  }

  template<class KTRAJ> void WireHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    ost << " WireHit constraining ";
    switch(wstate_.dimension_) {
      case WireHitState::none: default:
	ost << "Nothing";
	break;
      case WireHitState::time:
	ost << "Time";
	break;
      case WireHitState::distance:
	ost << "Distance";
	break;
      case WireHitState::both:
	ost << "Distance+Time";
	break;
    }
    ost << " LR Ambiguity " ;
    switch(wstate_.lrambig_) {
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
      if(active(WireHitState::time))
	ost << " Time Residual " << rresid_[WireHitState::time];
      if(active(WireHitState::distance))
	ost << " Distance Residual " << rresid_[WireHitState::distance];
      ost << std::endl;
    }
  }


} // KinKal namespace
#endif
