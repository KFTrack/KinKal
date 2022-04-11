#ifndef KinKal_StrawXing_hh
#define KinKal_StrawXing_hh
//
//  Describe the material effects of a kinematic trajectory crossing a straw
//  Used in the kinematic Kalman fit
//
#include "KinKal/Detector/ElementXing.hh"
#include "KinKal/Detector/StrawMaterial.hh"
#include "KinKal/Detector/StrawXingConfig.hh"
#include "KinKal/Detector/WireHit.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"

namespace KinKal {
  template <class KTRAJ> class StrawXing : public ElementXing<KTRAJ> {
    public:
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using EXING = ElementXing<KTRAJ>;
      using PCA = PiecewiseClosestApproach<KTRAJ,Line>;
      // construct from PCA
      StrawXing(PCA const& pca, StrawMaterial const& smat) : tpdata_(pca.tpData()), tprec_(pca.precision()), smat_(smat),
      sxconfig_(0.05*smat.strawRadius(),1.0),
      axis_(pca.sensorTraj()) {
        update(pca); }
      virtual ~StrawXing() {}
      // ElementXing interface
      void update(PKTRAJ const& pktraj) override;
      void update(PKTRAJ const& pktraj,MetaIterConfig const& miconfig) override;
      double time() const override { return tpdata_.particleToca(); }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // specific interface: this xing is based on PCA
      void update(PCA const& pca);
      // accessors
      ClosestApproachData const& closestApproach() const { return tpdata_; }
      double pcaPrecision() const { return tprec_; }
      StrawMaterial const& strawMaterial() const { return smat_; }
      StrawXingConfig const& config() const { return sxconfig_; }
    private:
      ClosestApproachData tpdata_; // result of most recent TPOCA
      double tprec_; // precision of TPOCA
      StrawMaterial const& smat_;
      StrawXingConfig sxconfig_;
      Line axis_; // straw axis, expressed as a timeline

      // should add state for displace wire TODO
  };

  template <class KTRAJ> void StrawXing<KTRAJ>::update(PCA const& pca) {
    if(pca.usable()){
      EXING::matXings().clear();
      tpdata_ = pca.tpData();
      smat_.findXings(tpdata_,sxconfig_,EXING::matXings());
    } else
      throw std::runtime_error("CA failure");
  }

  template <class KTRAJ> void StrawXing<KTRAJ>::update(PKTRAJ const& pktraj) {
    // use current xing time create a hint to the CA calculation: this speeds it up
    CAHint tphint(this->time(), this->time());
    PCA pca(pktraj,axis_,tphint,tprec_);
    update(pca);
  }

  template <class KTRAJ> void StrawXing<KTRAJ>::update(PKTRAJ const& pktraj,MetaIterConfig const& miconfig) {
    // search for an update to the xing configuration among this meta-iteration payload
    auto sxconfig = miconfig.findUpdater<StrawXingConfig>();
    if(sxconfig != 0) sxconfig_ = *sxconfig;
    update(pktraj);
  }

  template <class KTRAJ> void StrawXing<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost <<"Straw Xing time " << this->time();
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
