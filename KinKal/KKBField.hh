#ifndef KinKal_KKBField_hh
#define KinKal_KKBField_hh
//
// Class to correct for BField inhomogenity; adjust the parameters for the BField and momentum change
// The effect is located at a domain boundary
// This effect adds no information content or noise (presently), just transports the parameters 
//
#include "KinKal/KKEff.hh"
#include "KinKal/TDir.hh"
#include "KinKal/BField.hh"
#include "KinKal/BFieldUtils.hh"
#include <iostream>
#include <stdexcept>
#include <array>
#include <ostream>

namespace KinKal {
  template<class KTRAJ> class KKBField : public KKEff<KTRAJ> {
    public:
      using KKEFF = KKEff<KTRAJ>;
      using PKTRAJ = PKTraj<KTRAJ>;
      
      virtual double time() const override { return drange_.mid(); } // apply the correction at the middle of the range
      virtual bool isActive() const override { return active_ && bfcorr_ != KKConfig::nocorr; }
      virtual void update(PKTRAJ const& ref) override;
      virtual void update(PKTRAJ const& ref, MIConfig const& miconfig) override;
      virtual void print(std::ostream& ost=std::cout,int detail=0) const override;
      virtual void process(KKData& kkdata,TDir tdir) override;
      virtual void append(PKTRAJ& fit) override;
      DVEC const& effect() const { return dbeff_; }
      virtual ~KKBField(){}
      // create from the domain range, the effect, and the
      KKBField(BField const& bfield, PKTRAJ const& pktraj,TRange const& drange,KKConfig::BFieldCorr bfcorr) : 
	bfield_(bfield), drange_(drange), active_(false), bfcorr_(bfcorr) {} // not active until updated
	// accessors
	Vec3 deltaP() const { return Vec3(dp_[0], dp_[1], dp_[2]); } // translate to spatial vector
	TRange const& range() const { return drange_; }

    private:
      BField const& bfield_; // bfield
      SVec3 dp_; // change in momentum due to BField approximation
      TRange drange_; // extent of this effect.  The middle is at the transition point between 2 bfield domains (domain transition)
      DVEC dbint_; // integral effect of using bnom vs the full field over this effects range 
      PData dbeff_; // aggregate effect in parameter space of BField change, including BNom change
      bool active_; // activity state
      KKConfig::BFieldCorr bfcorr_; // type of correction to apply
  };

  template<class KTRAJ> void KKBField<KTRAJ>::process(KKData& kkdata,TDir tdir) {
    if(active_){
      // forwards; just append the effect's parameter change
      if(tdir == TDir::forwards) {
	kkdata.append(dbeff_);
      } else {
	// SUBTRACT the effect going backwards: covariance change is sign-independent
	PData reverse(dbeff_);
	reverse.parameters() *= -1.0;
      	kkdata.append(reverse);
      }
    }
    KKEffBase::setStatus(tdir,KKEffBase::processed);
  }

  template<class KTRAJ> void KKBField<KTRAJ>::update(PKTRAJ const& ref) {
    double etime = this->time();
    auto const& midtraj = ref.nearestPiece(etime);
    // compute parameter change due to integral of difference in BField vs BNom
    dbint_ = midtraj.dPardM(etime)*dp_;
    dbeff_.parameters() = dbint_;
    // add in the effect of changing BNom across this domain transition to parameters 
    if(bfcorr_ == KKConfig::variable){
      auto const& begtraj = ref.nearestPiece(drange_.begin());
      auto const& endtraj = ref.nearestPiece(drange_.end());
      dbeff_.parameters() += begtraj.dPardB(etime,endtraj.bnom()); // check sign FIXME!
    }
    // eventually include field map uncertainties in dbeff_ covariance TODO!
    KKEffBase::updateStatus();
  }

  template<class KTRAJ> void KKBField<KTRAJ>::update(PKTRAJ const& ref, MIConfig const& miconfig) {
    if(miconfig.updatebfcorr_){
      active_ = true;
      // integrate the fractional momentum change WRT this reference trajectory
      Vec3 dp =  BFieldUtils::integrate(bfield_, ref, drange_);
      dp_ = SVec3(dp.X(),dp.Y(),dp.Z()); //translate to SVec; this should be supported by SVector and GenVector
      //      std::cout << "Updating iteration " << miconfig.miter_ << " dP " << dp << std::endl;
    }
    update(ref);
  }

  template<class KTRAJ> void KKBField<KTRAJ>::append(PKTRAJ& fit) {
    if(active_){
      // make sure the piece is appendable
      if(fit.back().range().begin() > drange_.end()) throw std::invalid_argument("KKBField: Can't append piece");
      // adjust time if necessary
      double time = this->time()+ 1.0e-5; // slight buffer to make local piece selection more consistent
      double tlow = std::max(time,fit.back().range().begin() + 1.0e-5);
      TRange newrange(tlow,fit.range().end());
      KTRAJ newpiece(fit.back());
      newpiece.range() = newrange;
      // if we are using variable BField, update the parameters accordingly
      if(bfcorr_ == KKConfig::variable){
	Vec3 newbnom = bfield_.fieldVect(fit.position(drange_.end()));
	newpiece.setBNom(time,newbnom);
      }
      // adjust for the residual parameter change due to difference in bnom
      // don't double-count the effect due to bnom change; here we want just
      // the effect of the approximation of (piecewise) bnom vs the full field
      newpiece.params().parameters() += dbint_;
      fit.append(newpiece);
    }
  }

  template<class KTRAJ> void KKBField<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost << "KKBField " << static_cast<KKEff<KTRAJ>const&>(*this);
    ost << " dP " << dp_ << " effect " << dbeff_.parameters() << " domain range " << drange_ << std::endl;
  }

  template <class KTRAJ> std::ostream& operator <<(std::ostream& ost, KKBField<KTRAJ> const& kkmat) {
    kkmat.print(ost,0);
    return ost;
  }

}
#endif
