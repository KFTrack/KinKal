#ifndef KinKal_ShellXing_hh
#define KinKal_ShellXing_hh
//
//  Describe the effects of a kinematic trajectory crossing a thin shell of material defined by a surface
//  Used in the kinematic Kalman fit
//
#include "KinKal/Detector/ElementXing.hh"
#include "KinKal/Geometry/Surface.hh"
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
      void updateReference(KTRAJPTR const& ktrajptr) override;
      void updateState(MetaIterConfig const& config,bool first) override;
      Parameters params() const override;
      double time() const override { return tpca_.particleToca() + toff_; } // offset time WRT TOCA to avoid exact overlapp with the wire hit
      double transitTime() const override; // time to cross this element
      KTRAJ const& referenceTrajectory() const override { return tpca_.particleTraj(); }
      std::vector<MaterialXing>const&  matXings() const override { return mxings_; }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // accessors
    private:
      SURFPTR surf_; // surface
      MatEnv::DetMaterial const& mat_;
      Intersection inter_; // most recent intersection
      std::vector<MaterialXing> mxings_; // material xing
      double thick_; // shell thickness
      double varscale_; // variance scale, for annealing
      Parameters fparams_; // 1st-order parameter change for forwards time
  };

  template <class KTRAJ> ShellXing<KTRAJ>::ShellXing(SURFPTR surface, MatEnv::DetMaterial const& mat, Intersection inter, double thickness) :
    surf_(surface), mat_(mat), inter_(inter),thick_(thick),
    varscale_(1.0)
  {}

  template <class KTRAJ> void ShellXing<KTRAJ>::updateReference(KTRAJPTR const& ktrajptr) {
    // need a ptraj to re-compute the intersection, not clear what to do here FIXME
  }

  template <class KTRAJ> void ShellXing<KTRAJ>::updateState(MetaIterConfig const& miconfig,bool first) {
    if(first) {
      // search for an update to the xing configuration among this meta-iteration payload
      auto sxconfig = miconfig.findUpdater<ShellXingConfig>();
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
      double mommag = referenceTrajectory().momentum(time());
      // loop over the momentum change basis directions, adding up the effects on parameters from each
      for(int idir=0;idir<MomBasis::ndir; idir++) {
        auto mdir = static_cast<MomBasis::Direction>(idir);
        auto dir = referenceTrajectory().direction(time(),mdir);
        // project the momentum derivatives onto this direction
        DVEC pder = mommag*(dPdM*SVEC3(dir.X(), dir.Y(), dir.Z()));
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

  template <class KTRAJ> Parameters ShellXing<KTRAJ>::params() const {
    return fparams_;
  }

  template <class KTRAJ> double ShellXing<KTRAJ>::transitTime() const {
    return smat_.transitLength(tpca_.tpData())/tpca_.particleTraj().speed(tpca_.particleToca());
  }

  template <class KTRAJ> void ShellXing<KTRAJ>::print(std::ostream& ost,int detail) const {
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
