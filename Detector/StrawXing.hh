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
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      using EXING = ElementXing<KTRAJ>;
      using PCA = PiecewiseClosestApproach<KTRAJ,Line>;
      using CA = ClosestApproach<KTRAJ,Line>;
      // construct from PCA and material
      StrawXing(PCA const& pca, StrawMaterial const& smat);
      virtual ~StrawXing() {}
      // ElementXing interface
      void updateReference(KTRAJPTR const& ktrajptr) override;
      void updateState(MetaIterConfig const& config) override;
      double time() const override { return tpca_.particleToca() + toff_; } // offset time WRT TOCA to avoid exact overlapp with the wire hit
      KTRAJ const& referenceTrajectory() const override { return tpca_.particleTraj(); }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // accessors
      auto const& closestApproach() const { return tpca_; }
      auto const& strawMaterial() const { return smat_; }
      auto const& config() const { return sxconfig_; }
      auto precision() const { return tpca_.precision(); }
   private:
      Line axis_; // straw axis, expressed as a timeline
      StrawMaterial const& smat_;
      CA tpca_; // result of most recent TPOCA
      double toff_; // small time offset
      StrawXingConfig sxconfig_;
      // should add state for displaced wire/straw TODO
  };

  template <class KTRAJ> StrawXing<KTRAJ>::StrawXing(PCA const& pca, StrawMaterial const& smat) :
    axis_(pca.sensorTraj()),
    smat_(smat),
    tpca_(pca.localTraj(),axis_,pca.precision(),pca.tpData(),pca.dDdP(),pca.dTdP()),
    toff_(smat.strawRadius()*0.5/pca.particleTraj().speed(pca.particleToca())), // locate the effect to 1 side of the wire
    sxconfig_(0.05*smat.strawRadius(),1.0) // hardcoded values, should come from outside, FIXME
  {
//    this->updateReference(tpca_.particleTrajPtr());
  }

  template <class KTRAJ> void StrawXing<KTRAJ>::updateReference(KTRAJPTR const& ktrajptr) {
    CAHint tphint = tpca_.usable() ?  tpca_.hint() : CAHint(axis_.range().mid(),axis_.range().mid());
    tpca_ = CA(ktrajptr,axis_,tphint,precision());
    if(!tpca_.usable())throw std::runtime_error("WireHit TPOCA failure");
    // update the material effects
    smat_.findXings(tpca_.tpData(),sxconfig_,EXING::matXings());
 }

  template <class KTRAJ> void StrawXing<KTRAJ>::updateState(MetaIterConfig const& miconfig) {
    // search for an update to the xing configuration among this meta-iteration payload
    auto sxconfig = miconfig.findUpdater<StrawXingConfig>();
    if(sxconfig != 0){
      sxconfig_ = *sxconfig;
      smat_.findXings(tpca_.tpData(),sxconfig_,EXING::matXings());
    }
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
