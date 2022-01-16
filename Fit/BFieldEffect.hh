#ifndef KinKal_BFieldEffect_hh
#define KinKal_BFieldEffect_hh
//
// correct for the effect of BFieldMap inhomogenity; adjust the trajectory parameters using the BFieldMap
// This effect adds no information content or noise (presently), just transports the parameters
//
#include "KinKal/Fit/Effect.hh"
#include "KinKal/General/TimeDir.hh"
#include "KinKal/General/BFieldMap.hh"
#include "KinKal/Fit/Config.hh"
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
      bool active() const override { return bfcorr_; }
      void update(PKTRAJ const& ref) override;
      void update(PKTRAJ const& ref, MetaIterConfig const& miconfig) override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      void process(FitState& kkdata,TimeDir tdir) override;
      void append(PKTRAJ& fit) override;
      Parameters const& effect() const { return dbforw_; }
      virtual ~BFieldEffect(){}
      // disallow copy and equivalence
      BFieldEffect(BFieldEffect const& ) = delete;
      BFieldEffect& operator =(BFieldEffect const& ) = delete;
      // create from the domain range, the effect, and the
      BFieldEffect(Config const& config, BFieldMap const& bfield,TimeRange const& drange) :
        bfield_(bfield), drange_(drange), bfcorr_(config.bfcorr_) {}
      TimeRange const& range() const { return drange_; }

    private:
      BFieldMap const& bfield_; // bfield
      TimeRange drange_; // extent of this effect.  The middle is at the transition point between 2 bfield domains (domain transition)
      Parameters dbforw_; // aggregate effect in parameter space of BFieldMap change in the forwards direction
      bool bfcorr_; // apply correction or not
  };

  template<class KTRAJ> void BFieldEffect<KTRAJ>::process(FitState& kkdata,TimeDir tdir) {
    if(bfcorr_){
      // forwards; just append the effect's parameter change
      if(tdir == TimeDir::forwards) {
        kkdata.append(dbforw_);
      } else {
        // SUBTRACT the effect going backwards: covariance change is sign-independent
        Parameters reverse(dbforw_);
        reverse.parameters() *= -1.0;
        kkdata.append(reverse);
      }
    }
    KKEFF::setState(tdir,KKEFF::processed);
  }

  template<class KTRAJ> void BFieldEffect<KTRAJ>::update(PKTRAJ const& ref) {
    double etime = time();
    // compute parameter change due to changing BNom across this domain
    if(bfcorr_){
      auto const& begtraj = ref.nearestPiece(drange_.begin());
      auto const& endtraj = ref.nearestPiece(drange_.end());
      dbforw_.parameters() = begtraj.dPardB(etime,endtraj.bnom());
    }
    // eventually include field map uncertainties in dbforw_ covariance TODO!
    KKEFF::updateState();
  }

  template<class KTRAJ> void BFieldEffect<KTRAJ>::update(PKTRAJ const& ref, MetaIterConfig const& miconfig) {
    update(ref);
  }

  template<class KTRAJ> void BFieldEffect<KTRAJ>::append(PKTRAJ& fit) {
    if(bfcorr_){
      double etime = time();
      // make sure the piece is appendable
      if(fit.back().range().begin() > etime) throw std::invalid_argument("BFieldEffect: Can't append piece");
      TimeRange newrange(etime,std::max(fit.range().end(),drange_.end()));
      // copy the back piece of fit and set its range
      KTRAJ newpiece(fit.back());
      newpiece.range() = newrange;
      // update the parameters according to the change in bnom across this domain
      VEC3 newbnom = bfield_.fieldVect(fit.position3(drange_.end()));
      newpiece.setBNom(etime,newbnom);
      fit.append(newpiece);
    }
  }

  template<class KTRAJ> void BFieldEffect<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost << "BFieldEffect " << static_cast<Effect<KTRAJ>const&>(*this);
    ost << " effect " << dbforw_.parameters() << " domain range " << drange_ << std::endl;
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, BFieldEffect<KTRAJ> const& kkmat) {
    kkmat.print(ost,0);
    return ost;
  }

}
#endif
