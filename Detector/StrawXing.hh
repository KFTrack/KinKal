#ifndef KinKal_StrawXing_hh
#define KinKal_StrawXing_hh
//
//  Describe the material effects of a kinematic trajectory crossing a straw
//  Used in the kinematic Kalman fit
//
#include "KinKal/Detector/ElementXing.hh"
#include "KinKal/Detector/StrawMaterial.hh"
#include "KinKal/Detector/StrawXingConfig.hh"
#include "KinKal/Trajectory/SensorLine.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"

namespace KinKal {
  template <class KTRAJ> class StrawXing : public ElementXing<KTRAJ> {
    public:
      using PTRAJ = ParticleTrajectory<KTRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      using EXING = ElementXing<KTRAJ>;
      using PCA = PiecewiseClosestApproach<KTRAJ,SensorLine>;
      using CA = ClosestApproach<KTRAJ,SensorLine>;
      // construct from PCA and material
      StrawXing(PCA const& pca, StrawMaterial const& smat);
      virtual ~StrawXing() {}
      // ElementXing interface
      void updateReference(PTRAJ const& ptraj) override;
      void updateState(MetaIterConfig const& config,bool first) override;
      Parameters params() const override;
      double time() const override { return tpca_.particleToca() + toff_; } // offset time WRT TOCA to avoid exact overlapp with the wire hit
      double transitTime() const override; // time to cross this element
      KTRAJ const& referenceTrajectory() const override { return tpca_.particleTraj(); }
      std::vector<MaterialXing>const&  matXings() const override { return mxings_; }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // accessors
      auto const& closestApproach() const { return tpca_; }
      auto const& strawMaterial() const { return smat_; }
      auto const& config() const { return sxconfig_; }
      auto precision() const { return tpca_.precision(); }
    private:
      SensorLine axis_; // straw axis, expressed as a timeline
      StrawMaterial const& smat_;
      CA tpca_; // result of most recent TPOCA
      double toff_; // small time offset
      StrawXingConfig sxconfig_;
      double varscale_; // variance scale
      std::vector<MaterialXing> mxings_;
      Parameters fparams_; // parameter change for forwards time
  };

  template <class KTRAJ> StrawXing<KTRAJ>::StrawXing(PCA const& pca, StrawMaterial const& smat) :
    axis_(pca.sensorTraj()),
    smat_(smat),
    tpca_(pca.localTraj(),axis_,pca.precision(),pca.tpData(),pca.dDdP(),pca.dTdP()),
    toff_(smat.wireRadius()/pca.particleTraj().speed(pca.particleToca())), // locate the effect to 1 side of the wire to avoid overlap with hits
    varscale_(1.0)
  {}

  template <class KTRAJ> void StrawXing<KTRAJ>::updateReference(PTRAJ const& ptraj) {
    CAHint tphint = tpca_.usable() ?  tpca_.hint() : CAHint(axis_.timeAtMidpoint(),axis_.timeAtMidpoint());
    PCA pca(ptraj,axis_,tphint,precision());
    tpca_ = pca.localClosestApproach();
    if(!tpca_.usable())throw std::runtime_error("StrawXing TPOCA failure");
  }

  template <class KTRAJ> void StrawXing<KTRAJ>::updateState(MetaIterConfig const& miconfig,bool first) {
    if(first) {
      // search for an update to the xing configuration among this meta-iteration payload
      auto sxconfig = miconfig.findUpdater<StrawXingConfig>();
      if(sxconfig != 0){
        sxconfig_ = *sxconfig;
      }
      if(sxconfig_.scalevar_)
        varscale_ = miconfig.varianceScale();
      else
        varscale_ = 1.0;
    }
    smat_.findXings(tpca_.tpData(),sxconfig_,mxings_);
    // reset
    fparams_ = Parameters();
    if(mxings_.size() > 0){
      // compute the parameter effect for forwards time
      std::array<double,3> dmom = {0.0,0.0,0.0}, momvar = {0.0,0.0,0.0};
      this->materialEffects(dmom, momvar);
      // get the parameter derivative WRT momentum
      DPDV dPdM = referenceTrajectory().dPardM(time());
      // loop over the momentum change basis directions, adding up the effects on parameters from each
      for(int idir=0;idir<MomBasis::ndir; idir++) {
        auto mdir = static_cast<MomBasis::Direction>(idir);
        auto dir = referenceTrajectory().direction(time(),mdir);
        // project the momentum derivatives onto this direction
        DVEC pder = dPdM*SVEC3(dir.X(), dir.Y(), dir.Z());
        // convert derivative vector to a Nx1 matrix
        ROOT::Math::SMatrix<double,NParams(),1> dPdm;
        dPdm.Place_in_col(pder,0,0);
        // update the transport for this effect; Forward time propagation corresponds to energy loss
        fparams_.parameters() += pder*dmom[idir];
        // now the variance: this doesn't depend on time direction
        ROOT::Math::SMatrix<double, 1,1, ROOT::Math::MatRepSym<double,1>> MVar;
        MVar(0,0) = momvar[idir]*varscale_;
        fparams_.covariance() += ROOT::Math::Similarity(dPdm,MVar);
      }
    }
  }

  template <class KTRAJ> Parameters StrawXing<KTRAJ>::params() const {
    return fparams_;
  }

  template <class KTRAJ> double StrawXing<KTRAJ>::transitTime() const {
    return smat_.transitLength(tpca_.tpData())/tpca_.particleTraj().speed(tpca_.particleToca());
  }

  template <class KTRAJ> void StrawXing<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost <<"Straw Xing time " << this->time();
    if(detail > 0){
      for(auto const& mxing : mxings_){
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
