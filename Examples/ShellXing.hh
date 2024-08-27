#ifndef KinKal_ShellXing_hh
#define KinKal_ShellXing_hh
//
//  Describe the effects of a kinematic trajectory crossing a thin shell of material defined by a surface
//  Used in the kinematic Kalman fit
//
#include "KinKal/Detector/ElementXing.hh"
#include "KinKal/Geometry/Surface.hh"
#include "KinKal/Geometry/ParticleTrajectoryIntersect.hh"
#include "Offline/Mu2eKinKal/inc/ShellXingUpdater.hh"
#include "KinKal/MatEnv/DetMaterial.hh"

namespace KinKal {
  template <class KTRAJ,class SURF> class ShellXing : public ElementXing<KTRAJ> {
    public:
      using PTRAJ = ParticleTrajectory<KTRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      using EXING = ElementXing<KTRAJ>;
      using PCA = PiecewiseClosestApproach<KTRAJ,SensorLine>;
      using CA = ClosestApproach<KTRAJ,SensorLine>;
      using SURFPTR = std::shared_ptr<SURF>;
      // construct from a surface, material, intersection, and transverse thickness
      ShellXing(SURFPTR surface, MatEnv::DetMaterial const& mat, Intersection inter, double thickness)
      virtual ~ShellXing() {}
      // ElementXing interface
      void updateReference(PTRAJ const& ptraj) override;
      void updateState(MetaIterConfig const& config,bool first) override;
      Parameters params() const override;
      double time() const override { return inter_.time(); }
      KTRAJ const& referenceTrajectory() const override { return *reftrajptr_; }
      std::vector<MaterialXing>const&  matXings() const override { return mxings_; }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // accessors
    private:
      SURFPTR surf_; // surface
      MatEnv::DetMaterial const& mat_;
      Intersection inter_; // most recent intersection
      KTRAJPTR reftrajptr_; // reference trajectory
      std::vector<MaterialXing> mxings_; // material xing
      double thick_; // shell thickness
      double tol_; // tolerance for intersection
      Parameters fparams_; // 1st-order parameter change for forwards time
      ShellXingUpdater sxconfig_; // note this must come from an updater during processing
      double varscale_; // cache
  };

  template <class KTRAJ> ShellXing<KTRAJ>::ShellXing(SURFPTR surface, MatEnv::DetMaterial const& mat, Intersection inter,
      KTTRAJPTR reftrajptr, double thickness, double tol) :
    surf_(surface), mat_(mat), inter_(inter), reftrajptr_(reftrajptr), thick_(thick),tol_(tol), varscale_(1.0)
  {}

  template <class KTRAJ> void ShellXing<KTRAJ>::updateReference(PTRAJ const& ptraj) {
    // re-intersect with the surface, taking the current time as start and range from the current piece (symmetrized)
    double delta = 0.5*reftraj->range().range(); // this factor may need tuning: too small and it may miss the intersection, too large and it may jump to a different intersection TODO
    TimeRange irange(inter_.time_-delta, inter_.time_+delta);
    inter_ = intersect(ptraj, &surf_, irange,tol_);
    reftrajptr_ = ptraj.nearestTraj(inter_.time_); // I may need to protect against the case the time is crazy (failed intersection) TODO
  }

  template <class KTRAJ> void ShellXing<KTRAJ>::updateState(MetaIterConfig const& miconfig,bool first) {
    if(first) {
      // search for an update to the xing configuration among this meta-iteration payload
      auto sxconfig = miconfig.findUpdater<ShellXingUpdater>();
      if(sxconfig != 0) sxconfig_ = *sxconfig;
      if(sxconfig_.scalevar_) varscale_ = miconfig.varianceScale();
      // find the material xings from gas, straw wall, and wire
      auto cad = ca_.tpData();
      if(shptr_ && shptr_->hitState().active()){
        // if we have an associated hit, overwrite the DOCA and DOCAVAR using the drift info, which is much more accurate
        auto dinfo = shptr_->fillDriftInfo();
        cad.doca_ = dinfo.rDrift_;
        cad.docavar_ = dinfo.unsignedDriftVar();
      }
      mxings_.clear();
      // check if we are on the surface
      if(inter_.onsurface_ && inter_.inbounds_){
        // compute the material
      }
    }
    if(mxings_.size() > 0){
      fparams_ = this->parameterChange(varscale_);
    } else {
      // reset
      fparams_ = Parameters();
    }
  }

  template <class KTRAJ> Parameters ShellXing<KTRAJ>::params() const {
    return fparams_;
  }

  template <class KTRAJ> void ShellXing<KTRAJ>::print(std::ostream& ost,int detail) const {
    ost <<"Shell Xing time " << this->time();
    ost << std::endl;
  }

}
#endif
