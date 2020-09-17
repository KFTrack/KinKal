#ifndef KinKal_BFieldEffect_hh
#define KinKal_BFieldEffect_hh
//
// correct for the effect of BFieldMap inhomogenity; adjust the trajectory parameters using the BFieldMap
// This effect adds no information content or noise (presently), just transports the parameters 
//
#include "KinKal/Effect.hh"
#include "KinKal/TimeDir.hh"
#include "KinKal/BFieldMap.hh"
#include "KinKal/BFieldUtils.hh"
#include "KinKal/Config.hh"
#include <iostream>
#include <stdexcept>
#include <array>
#include <ostream>

namespace KinKal {
  template<class KTRAJ> class BFieldEffect : public Effect<KTRAJ> {
    public:
      using KKEFF = Effect<KTRAJ>;
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      
      double time() const override { return drange_.mid(); } // apply the correction at the middle of the range
      bool isActive() const override { return active_ && bfcorr_ != Config::nocorr; }
      void update(PKTRAJ const& ref) override;
      void update(PKTRAJ const& ref, MetaIterConfig const& miconfig) override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      void process(FitState& kkdata,TimeDir tdir) override;
      void append(PKTRAJ& fit) override;
      Parameters const& effect() const { return dbeff_; }
      virtual ~BFieldEffect(){}
      // disallow copy and equivalence
      BFieldEffect(BFieldEffect const& ) = delete; 
      BFieldEffect& operator =(BFieldEffect const& ) = delete; 
      // create from the domain range, the effect, and the
      BFieldEffect(Config const& config, PKTRAJ const& pktraj,TimeRange const& drange) : 
	bfield_(config.bfield_), drange_(drange), active_(false), bfcorr_(config.bfcorr_) {}
      VEC3 deltaP() const { return VEC3(dp_[0], dp_[1], dp_[2]); } // translate to spatial vector
      TimeRange const& range() const { return drange_; }

    private:
      BFieldMap const& bfield_; // bfield
      SVEC3 dp_; // change in momentum due to BFieldMap approximation
      TimeRange drange_; // extent of this effect.  The middle is at the transition point between 2 bfield domains (domain transition)
      DVEC dbint_; // integral effect of using bnom vs the full field over this effects range 
      Parameters dbeff_; // aggregate effect in parameter space of BFieldMap change, including BNom change
      bool active_; // activity state
      Config::BFCorr bfcorr_; // type of BFieldMap map correction to apply
  };

  template<class KTRAJ> void BFieldEffect<KTRAJ>::process(FitState& kkdata,TimeDir tdir) {
    if(active_){
      // forwards; just append the effect's parameter change
      if(tdir == TimeDir::forwards) {
	kkdata.append(dbeff_);
      } else {
	// SUBTRACT the effect going backwards: covariance change is sign-independent
	Parameters reverse(dbeff_);
	reverse.parameters() *= -1.0;
      	kkdata.append(reverse);
      }
    }
    KKEFF::setStatus(tdir,KKEFF::processed);
  }

  template<class KTRAJ> void BFieldEffect<KTRAJ>::update(PKTRAJ const& ref) {
    double etime = this->time();
    auto const& midtraj = ref.nearestPiece(etime);
    // compute parameter change due to integral of difference in BFieldMap vs BNom
    if(bfcorr_ == Config::fixed || bfcorr_ == Config::both){
      dbint_ = midtraj.dPardM(etime)*dp_;
    } else {
      dbint_ = DVEC();
    }
    dbeff_.parameters() = dbint_;
    // add in the effect of changing BNom across this domain transition to parameters 
    if(bfcorr_ == Config::variable || bfcorr_ == Config::both){
      auto const& begtraj = ref.nearestPiece(drange_.begin());
      auto const& endtraj = ref.nearestPiece(drange_.end());
      dbeff_.parameters() += begtraj.dPardB(etime,endtraj.bnom()); // check sign FIXME!
    }
    // eventually include field map uncertainties in dbeff_ covariance TODO!
    KKEFF::updateStatus();
  }

  template<class KTRAJ> void BFieldEffect<KTRAJ>::update(PKTRAJ const& ref, MetaIterConfig const& miconfig) {
    if(miconfig.updatebfcorr_){
      active_ = true;
      // integrate the fractional momentum change WRT this reference trajectory
      VEC3 dp =  BFieldUtils::integrate(bfield_, ref, drange_);
      dp_ = SVEC3(dp.X(),dp.Y(),dp.Z()); //translate to SVec; this should be supported by SVector and GenVector
      //      std::cout << "Updating iteration " << miconfig.miter_ << " dP " << dp << std::endl;
    } else {
      active_ = false;
    }
    update(ref);
  }

  template<class KTRAJ> void BFieldEffect<KTRAJ>::append(PKTRAJ& fit) {
    if(active_){
      // make sure the piece is appendable
      if(fit.back().range().begin() > drange_.end()) throw std::invalid_argument("BFieldEffect: Can't append piece");
      // adjust time if necessary
      double time = this->time()+ 1.0e-5; // slight buffer to make local piece selection more consistent
      double tlow = std::max(time,fit.back().range().begin() + 1.0e-5);
      TimeRange newrange(tlow,fit.range().end());
      KTRAJ newpiece(fit.back());
      newpiece.range() = newrange;
      // if we are using variable BFieldMap, update the parameters accordingly
      if(bfcorr_ == Config::variable || bfcorr_ == Config::both){
	VEC3 newbnom = bfield_.fieldVect(fit.position3(drange_.end()));
	newpiece.setBNom(time,newbnom);
      }
      // adjust for the residual parameter change due to difference in bnom
      // don't double-count the effect due to bnom change; here we want just
      // the effect of the approximation of (piecewise) bnom vs the full field
      newpiece.params().parameters() += dbint_;
      fit.append(newpiece);
    }
  }

  template<class KTRAJ> void BFieldEffect<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost << "BFieldEffect " << static_cast<Effect<KTRAJ>const&>(*this);
    ost << " dP " << dp_ << " effect " << dbeff_.parameters() << " domain range " << drange_ << std::endl;
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, BFieldEffect<KTRAJ> const& kkmat) {
    kkmat.print(ost,0);
    return ost;
  }

}
#endif
