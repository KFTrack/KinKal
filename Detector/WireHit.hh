#ifndef KinKal_WireHit_hh
#define KinKal_WireHit_hh
//
//  class representing a drift wire measurement.  Implemented using PTCA between the particle traj and the wire
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/Detector/Hit.hh"
#include "KinKal/Trajectory/Residual.hh"
#include "KinKal/Detector/WireCell.hh"
#include "KinKal/Detector/ElementXing.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/General/LRAmbig.hh"
#include "KinKal/Detector/BFieldMap.hh"
#include <stdexcept>
namespace KinKal {
// struct for updating wire hits; this is just parameters, but could be methods as well
  struct WireHitUpdater {
    double mindoca_; // minimum DOCA value to set an ambiguity
    double maxdoca_; // maximum DOCA to still use a hit
    WireHitUpdater(double mindoca,double maxdoca) : mindoca_(mindoca), maxdoca_(maxdoca) {}
  };

  template <class KTRAJ> class WireHit : public Hit<KTRAJ> {
    public:
      using HIT = Hit<KTRAJ>;
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using PTCA = PiecewiseClosestApproach<KTRAJ,Line>;
      using EXING = ElementXing<KTRAJ>;
      using EXINGPTR = std::shared_ptr<EXING>;
      
      // Hit interface overrrides
      Weights weight() const override;
      double chi(Parameters const& pdata) const override;
      double time() const override { return rresid_.time(); }
      void update(PKTRAJ const& pktraj, MetaIterConfig const& config) override;
      void update(PKTRAJ const& pktraj) override;
      bool isActive() const override { return active_; }
      EXINGPTR const& detXingPtr() const override { return dxing_; }
      unsigned nDOF() const override { return active_ ? 1 : 0; }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      
      Line const& wire() const { return wire_; }
      // set the null variance given the min DOCA used to assign LR ambiguity.  This assumes a flat DOCA distribution
      void setNullVar(double mindoca) { nullvar_ = mindoca*mindoca/12.0; }
      void setAmbig(LRAmbig newambig) { ambig_ = newambig; }
      WireHit(BFieldMap const& bfield, Line const& wire, WireCell const& cell, EXINGPTR const& dxing, LRAmbig ambig=LRAmbig::null);
      virtual ~WireHit(){}
      LRAmbig ambig() const { return ambig_; }
      WireCell const& cell() const { return cell_; }
      Residual const& refResidual() const { return rresid_; }
      Parameters const& refParams() const { return rparams_; }
    private:
      void setResidual(PTCA const& tpoca);
      Line wire_; // local linear approximation to the wire of this hit.  The range describes the active wire length
      WireCell const& cell_; // cell description
      double nullvar_; // variance of the error in space for null ambiguity
      LRAmbig ambig_; // current ambiguity assignment: can change during a fit
      bool active_; // active or not (pat. rec. tool)
      BFieldMap const& bfield_;
      EXINGPTR dxing_; // material xing
      // caches used in processing
      Residual rresid_; // residual WRT most recent reference parameters
      Parameters rparams_; // reference parameters
      double precision_; // current precision
  };

  template <class KTRAJ> WireHit<KTRAJ>::WireHit(BFieldMap const& bfield, Line const& wire, WireCell const& cell, EXINGPTR const& dxing, LRAmbig ambig) : 
    wire_(wire), cell_(cell), ambig_(ambig), active_(true), bfield_(bfield), dxing_(dxing), precision_(1e-6) { setNullVar(cell_.size()); }

  template <class KTRAJ> Weights WireHit<KTRAJ>::weight() const {
    if(isActive()){
      // convert derivatives to a Nx1 matrix (for root)
      ROOT::Math::SMatrix<double,NParams(),1> dRdPM;
      dRdPM.Place_in_col(rresid_.dRdP(),0,0);
      // convert the variance into a 1X1 matrix
      ROOT::Math::SMatrix<double, 1,1, ROOT::Math::MatRepSym<double,1>> RVarM;
      // weight by inverse variance
      double tvar = rresid_.variance();
      RVarM(0,0) = 1.0/tvar;
      // expand these into the weight matrix
      DMAT wmat = ROOT::Math::Similarity(dRdPM,RVarM);
      // translate residual value into weight vector WRT the reference parameters
      // sign convention reflects resid = measurement - prediction
      DVEC wvec = wmat*rparams_.parameters() + rresid_.dRdP()*rresid_.value()/tvar;
      return Weights(wvec,wmat);
    } else
      return Weights();
  }

  template <class KTRAJ> void WireHit<KTRAJ>::update(PKTRAJ const& pktraj) {
    // compute PTCA.  Use the wire middle as hint
    CAHint tphint(wire_.range().mid(),wire_.range().mid());
    if(rresid_.tPoca().usable())
      tphint = CAHint(rresid_.tPoca().particleToca(),rresid_.tPoca().sensorToca());
    PTCA tpoca(pktraj,wire_,tphint,precision_);
    setResidual(tpoca);
    rparams_ = pktraj.nearestPiece(rresid_.time()).params();
  }

  template <class KTRAJ> void WireHit<KTRAJ>::update(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) {
    precision_ = miconfig.tprec_;
    // update to move to the new trajectory
    update(pktraj);
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
      double doca = rresid_.tPoca().doca();
      if(fabs(doca) > whupdater->mindoca_){
	LRAmbig newambig = doca > 0.0 ? LRAmbig::right : LRAmbig::left;
	setAmbig(newambig);
      } else {
	setAmbig(LRAmbig::null);
	setNullVar(std::min(cell_.size(),whupdater->mindoca_));
      }
      // decide if the hit is consistent with this track, and if not disable/enable it.
      // for now, just look at DOCA, but could look at other information.  This should be
      // delegated to the hit updater as a function FIXME!
      active_ = fabs(doca) < whupdater->maxdoca_;
    // now update again in case the hit changed
      update(pktraj);
    }
  // OK if no updater is found, hits may be frozen this meta-iteration
  }

  template <class KTRAJ> void WireHit<KTRAJ>::setResidual(PTCA const& tpoca) {
    if(tpoca.usable()){
      // translate PTCA to residual
      VEC3 bvec = bfield_.fieldVect(tpoca.particlePoca().Vect());
      auto pdir = bvec.Cross(wire_.dir()).Unit(); // direction perp to wire and BFieldMap
      VEC3 dvec = tpoca.delta().Vect();
      double phi = asin(double(dvec.Unit().Dot(pdir)));
      // must use absolute DOCA to call distanceToTime
      POL2 drift(fabs(tpoca.doca()), phi);
      double tdrift, tdvar, vdrift;
      cell_.distanceToTime(drift, tdrift, tdvar, vdrift);
      // Use ambiguity to convert drift time to a time difference.   null ambiguity means ignore drift time
      auto iambig = static_cast<std::underlying_type<LRAmbig>::type>(ambig_);
      double dsign = copysign(1.0,iambig*tpoca.lSign()); // overall sign is the product of ambiguity and doca sign
      double dt = tpoca.deltaT()-tdrift*dsign;
      if(ambig_ != LRAmbig::null){ 
	// residual is in time, so unit dependendence on time, distance dependence is the local drift velocity
	DVEC dRdP = tpoca.dDdP()*iambig/vdrift - tpoca.dTdP(); 
	rresid_ = Residual(Residual::dtime,tpoca.tpData(),dt,tdvar,dRdP);
      } else {
	// interpret DOCA against the wire directly as the residual.  
	rresid_ = Residual(Residual::distance,tpoca.tpData(),-tpoca.doca(),nullvar_,tpoca.dDdP());
      }
    } else
      throw std::runtime_error("TCA failure");
  }

  template <class KTRAJ> double WireHit<KTRAJ>::chi(Parameters const& pdata) const {
    double retval(0.0);
    if(isActive()) {
      // compute the difference between these parameters and the reference parameters
      DVEC dpvec = pdata.parameters() - rparams_.parameters(); 
      // use the differnce to 'correct' the reference residual to be WRT these parameters
      double uresid = rresid_.value() - ROOT::Math::Dot(dpvec,rresid_.dRdP());
      // project the parameter covariance into a residual space variance
      double rvar = ROOT::Math::Similarity(rresid_.dRdP(),pdata.covariance());
      // add the measurement variance
      rvar +=  rresid_.variance();
      // chi is the ratio of these
      retval = uresid/sqrt(rvar);
    }
    return retval;
  }

  template<class KTRAJ> void WireHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    if(this->isActive())
      ost<<"Active ";
    else
      ost<<"Inactive ";
    ost << " WireHit LRAmbig " << this-> ambig() << " " << std::endl;
  }

}
#endif
