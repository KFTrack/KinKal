#ifndef KinKal_StrawXing_hh
#define KinKal_StrawXing_hh
//
//  Describe the material effects of a kinematic trajectory crossing a straw
//  Used in the kinematic Kalman fit
//
#include "KinKal/Detector/ElementXing.hh"
#include "KinKal/Detector/StrawMaterial.hh"
#include "KinKal/Detector/WireHit.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"

namespace KinKal {
  template <class KTRAJ> class StrawXing : public ElementXing<KTRAJ> {
    public:
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using EXING = ElementXing<KTRAJ>;
      using PTCA = PiecewiseClosestApproach<KTRAJ,Line>;
      using STRAWHIT = WireHit<KTRAJ>;
      using STRAWHITPTR = std::shared_ptr<STRAWHIT>;
      // construct from PTCA
      StrawXing(PTCA const& tpoca, StrawMaterial const& smat) : EXING(tpoca.particleToca()) , tprec_(tpoca.precision()), smat_(smat),
      sxconfig_(0.05*smat.strawRadius(),1.0),
      axis_(tpoca.sensorTraj()) {
	update(tpoca); }
      virtual ~StrawXing() {}
      // ElementXing interface
      void update(PKTRAJ const& pktraj,MetaIterConfig const& miconfig) override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
     // specific interface: this xing is based on PTCA
      void update(PTCA const& tpoca);
      // accessors
      StrawMaterial const& strawMaterial() const { return smat_; }
      StrawXingConfig const& config() const { return sxconfig_; }
    private:
      double tprec_;
      StrawMaterial const& smat_;
      StrawXingConfig sxconfig_;
      Line axis_; // straw axis, expressed as a timeline
      // should add state for displace wire TODO
  };

  template <class KTRAJ> void StrawXing<KTRAJ>::update(PTCA const& tpoca) {
    if(tpoca.usable()){
      EXING::matXings().clear();
      smat_.findXings(tpoca.tpData(),sxconfig_,EXING::matXings());
      EXING::crossingTime() = tpoca.particleToca();
    } else
      throw std::runtime_error("CA failure");
  }

  template <class KTRAJ> void StrawXing<KTRAJ>::update(PKTRAJ const& pktraj,MetaIterConfig const& miconfig) {
  // search for an update to the xing configuration among this meta-iteration payload
    const StrawXingConfig* sxconfig(0);
    for(auto const& uparams : miconfig.updaters_){
      auto const* sxc = std::any_cast<StrawXingConfig>(&uparams);
      if(sxc != 0){
	if(sxconfig !=0) throw std::invalid_argument("Multiple SimpleWireHitUpdaters found");
	sxconfig = sxc;
      }
    }
    if(sxconfig != 0) sxconfig_ = *sxconfig;
    // use current xing time create a hint to the CA calculation: this speeds it up
    CAHint tphint(EXING::crossingTime(), EXING::crossingTime());
    PTCA tpoca(pktraj,axis_,tphint,tprec_);
    update(tpoca);
  }

  template <class KTRAJ> void StrawXing<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost <<"Straw Xing time " << this->crossingTime();
    if(detail > 0){
      for(auto const& mxing : this->matXings()){
	ost << " " << mxing.dmat_.name() << " pathLen " << mxing.plen_;
      }
    }
    if(detail > 1){
      ost << " Axis ";
      axis_.print(ost,0);
    }
    ost << std::endl;
  }

}
#endif
