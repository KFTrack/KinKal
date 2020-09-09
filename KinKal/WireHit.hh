#ifndef KinKal_WireHit_hh
#define KinKal_WireHit_hh
//
//  class representing a drift wire measurement.  Implemented using PTPOCA between the particle traj and the wire
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/DetectorHit.hh"
#include "KinKal/DistanceToTime.hh"
#include "KinKal/DetectorXing.hh"
#include "KinKal/Line.hh"
#include "KinKal/PieceClosestApproach.hh"
#include "KinKal/LRAmbig.hh"
#include "KinKal/BFieldMap.hh"
#include <stdexcept>
namespace KinKal {
// struct for updating wire hits; this is just parameters, but could be methods as well
  struct WireHitUpdater {
    double mindoca_; // minimum DOCA value to set an ambiguity
    double maxdoca_; // maximum DOCA to still use a hit
    WireHitUpdater(double mindoca,double maxdoca) : mindoca_(mindoca), maxdoca_(maxdoca) {}
  };

  template <class KTRAJ> class WireHit : public DetectorHit<KTRAJ> {
    public:
      using THIT = DetectorHit<KTRAJ>;
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using PTPOCA = PieceClosestApproach<KTRAJ,Line>;
      using DXING = DetectorXing<KTRAJ>;
      using DXINGPTR = std::shared_ptr<DXING>;
      
      // Hit interface overrrides
      virtual void resid(PKTRAJ const& pktraj, Residual& resid, double precision) const override;
      void resid(PTPOCA const& tpoca, Residual& resid) const; // actual implementation of resid uses PTPOCA
      virtual void update(PKTRAJ const& pktraj, MetaIterConfig const& config, Residual& resid) override;
      virtual unsigned nDOF() const override { return 1; }
      double cellSize() const { return csize_; } // approximate transverse cell size, used to set null variance
// construct from a DistanceToTime relationship; BFieldMap is needed to compute ExB effects
      Line const& wire() const { return wire_; }
      // set the null variance given the min DOCA used to assign LR ambiguity.  This assumes a flat DOCA distribution
      void setNullVar(double mindoca) { nullvar_ = mindoca*mindoca/3.0; }
      void setAmbig(LRAmbig newambig) { ambig_ = newambig; }
      WireHit(DXINGPTR const& dxing, BFieldMap const& bfield, Line const& wire, DistanceToTime const& d2t, double csize,LRAmbig ambig=LRAmbig::null) : 
	THIT(dxing,true), wire_(wire), d2t_(d2t), csize_(csize), ambig_(ambig), bfield_(bfield) { setNullVar(csize_); }
      virtual ~WireHit(){}
      LRAmbig ambig() const { return ambig_; }
      DistanceToTime const& d2T() const { return d2t_; }
    private:
      Line wire_; // local linear approximation to the wire of this hit.  The range describes the active wire length
      DistanceToTime const& d2t_; // distance to time relationship for drift in this cell
      double csize_; // transverse cell size in mm
      double nullvar_; // variance of the error in space for null ambiguity
      LRAmbig ambig_; // current ambiguity assignment: can change during a fit
      BFieldMap const& bfield_;
  };

  template <class KTRAJ> void WireHit<KTRAJ>::resid(PKTRAJ const& pktraj, Residual& residual,double precision) const {
    // compute PTPOCA.  Use the wire middle as hint
    CAHint tphint(wire_.range().mid(),wire_.range().mid());
    PTPOCA tpoca(pktraj,wire_,tphint,precision);
    resid(tpoca,residual);
  }

  template <class KTRAJ> void WireHit<KTRAJ>::update(PKTRAJ const& pktraj, MetaIterConfig const& miconfig, Residual& residual ) {
    // find PTPOCA, using previous residual as hint, or the wire itself if not
    CAHint tphint(wire_.range().mid(),wire_.range().mid());
    if(residual.tPoca().usable())
      tphint = CAHint(residual.tPoca().particleToca(),residual.tPoca().sensorToca());
    PTPOCA tpoca(pktraj,wire(),tphint,miconfig.tprec_);
    // find the wire hit updater in the update params.  If there are more than 1 throw 
    const WireHitUpdater* whupdater(0);
    for(auto const& uparams : miconfig.hitupdaters_){
      auto const* whu = std::any_cast<WireHitUpdater>(&uparams);
      if(whu != 0){
	if(whupdater !=0) throw std::invalid_argument("Multiple WireHitUpdaters found");
	whupdater = whu;
      }
    }
    if(whupdater != 0){
        // use DOCA to set the ambiguity
      if(fabs(tpoca.doca()) > whupdater->mindoca_){
	LRAmbig newambig = tpoca.doca() > 0.0 ? LRAmbig::right : LRAmbig::left;
	setAmbig(newambig);
      } else {
	setAmbig(LRAmbig::null);
	setNullVar(std::min(cellSize(),whupdater->mindoca_));
      }
      // decide if the hit is consistent with this track, and if not disable/enable it.
      // for now, just look at DOCA, but could use tension too TODO!
      THIT::setActivity(fabs(tpoca.doca()) < whupdater->maxdoca_);
    } // allow no updater: hits may be frozen this meta-iteration
    // compute the residual
    resid(tpoca,residual);
  }

  template <class KTRAJ> void WireHit<KTRAJ>::resid(PTPOCA const& tpoca, Residual& resid) const {
    if(tpoca.usable()){
      // translate PTPOCA to residual
      if(ambig_ != LRAmbig::null){ 
	auto iambig = static_cast<std::underlying_type<LRAmbig>::type>(ambig_);
	// convert DOCA to wire-local polar coordinates.  This defines azimuth WRT the B field for ExB effects
	double rho = tpoca.doca()*iambig; // this is allowed to go negative
	VEC3 bvec = bfield_.fieldVect(tpoca.particlePoca().Vect());
	auto pdir = bvec.Cross(wire_.dir()).Unit(); // direction perp to wire and BFieldMap
	VEC3 dvec = tpoca.delta().Vect();
	double phi = asin(double(dvec.Unit().Dot(pdir)));
	POL2 drift(rho, phi);
	double tdrift, tdvar, vdrift;
	d2T().distanceToTime(drift, tdrift, tdvar, vdrift);
	// residual is in time, so unit dependendence on time, distance dependence is the local drift velocity
	DVEC dRdP = tpoca.dDdP()*iambig/vdrift - tpoca.dTdP(); 
	resid = Residual(Residual::dtime,tpoca.tpData(),tpoca.deltaT()-tdrift,tdvar,dRdP);
      } else {
	// interpret DOCA against the wire directly as the residual.  There is no direct time dependence in this case
	// residual is in space, so unit dependendence on distance, none on time
	resid = Residual(Residual::distance,tpoca.tpData(),-tpoca.doca(),nullvar_,tpoca.dDdP());
      }
    } else
      throw std::runtime_error("TPOCA failure");
  }

}
#endif
