#ifndef KinKal_ElementXing_hh
#define KinKal_ElementXing_hh
//
//  Describe the material effects of a particle crossing a detector element (piece of the detector)
//  Used in the kinematic Kalman fit
//
#include "KinKal/General/MomBasis.hh"
#include "KinKal/Detector/MaterialXing.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
#include "KinKal/Fit/MetaIterConfig.hh"
#include <vector>
#include <array>
#include <ostream>

namespace KinKal {
  template <class KTRAJ> class ElementXing {
    public:
      using PTRAJ = ParticleTrajectory<KTRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      ElementXing() {}
      virtual ~ElementXing() {}
      virtual void updateReference(PTRAJ const& ptraj) = 0; // update the trajectory reference
      virtual void updateState(MetaIterConfig const& config,bool first) =0; // update the state according to this meta-config
      virtual Parameters params() const =0; // parameter change induced by this element crossing WRT the reference parameters going forwards in time
      virtual double time() const=0; // time the particle crosses thie element
      virtual double transitTime() const=0; // time to cross this element
      virtual KTRAJ const& referenceTrajectory() const =0; // trajectory WRT which the xing is defined
      virtual std::vector<MaterialXing>const&  matXings() const =0; // Effect of each physical material component of this detector element on this trajectory
      virtual void print(std::ostream& ost=std::cout,int detail=0) const =0;
      virtual bool active() const =0;
      // momentum change and variance increase associated with crossing this element forwards in time, in spatial basis
      void momentumChange(SVEC3& dmom, SMAT& dmomvar) const;
      // parameter change associated with crossing this element forwards in time
      Parameters parameterChange(double varscale=1.0) const;
      // cumulative material effect from all the materials in this element crossing going forwards in time, expressed as scalars
      void materialEffects(double& dmom, double& paramomvar, double& perpmomvar) const;
      // sum radiation fraction
      double radiationFraction() const;
    private:
  };

  template <class KTRAJ> void ElementXing<KTRAJ>::momentumChange(SVEC3& dmom, SMAT& dmomvar) const {
    // compute the parameter effect for forwards time
    double dm, paramomvar, perpmomvar;
    materialEffects(dm, paramomvar,perpmomvar);
    // momentum change in forward time direction /mdue to energy loss; this is along the momentum
    auto momdir = referenceTrajectory().direction(time());
    dmom = dm*SVEC3(momdir.X(),momdir.Y(),momdir.Z());
    // now update the covariance; this includes smearing from energy straggling and multiple scattering
    // move variances into a matrix
    ROOT::Math::SVector<double, 6>  varvec(paramomvar, 0, perpmomvar, 0, 0, perpmomvar);
    SMAT mmvar(varvec);
    // loop over the momentum change basis directions and create the transform matrix between that and global Cartesian basis
    SSMAT dmdxyz; // momentum basis -> Cartesian conversion matrix
    for(int idir=0;idir<MomBasis::ndir; idir++) {
      auto mdir = referenceTrajectory().direction(time(),static_cast<MomBasis::Direction>(idir));
      SVEC3 vmdir(mdir.X(), mdir.Y(), mdir.Z());
      dmdxyz.Place_in_col(vmdir,0,idir);
    }
    // return variance in global Cartesian coordinates
    dmomvar = ROOT::Math::Similarity(dmdxyz,mmvar);
  }

  template <class KTRAJ> Parameters ElementXing<KTRAJ>::parameterChange(double varscale) const {
    // compute this xing's effect on momentum in global Cartesian
    SVEC3 dmom;
    SMAT dmomvar;
    momentumChange(dmom,dmomvar);
    // convert that to parameter space
    DPDV dPdM = referenceTrajectory().dPardM(time());
    auto dmomp = dPdM*dmom;
    // scale covariance as needed
    dmomvar*= varscale;
    auto dmompvar = ROOT::Math::Similarity(dPdM,dmomvar);
    return Parameters(dmomp,dmompvar);
  }

  template <class KTRAJ> void ElementXing<KTRAJ>::materialEffects(double& dmom, double& paramomvar, double& perpmomvar) const {
    // accumulate the change in energy and scattering angle variance from the material components
    double E = referenceTrajectory().energy();
    double mom = referenceTrajectory().momentum();
    double mass = referenceTrajectory().mass();
    // loop over individual materials and accumulate their effects
    double dE(0.0), dEVar(0.0), scatvar(0.0);
    for(auto const& mxing : matXings()){
      dE += mxing.dmat_.energyLoss(mom,mxing.plen_,mass);
      dEVar += mxing.dmat_.energyLossVar(mom,mxing.plen_,mass);
      scatvar += mxing.dmat_.scatterAngleVar(mom,mxing.plen_,mass);
    }
    // convert energy change to change in momentum
    double dmdE = E/mom;
    dmom = dE*dmdE;
    paramomvar = dEVar*dmdE*dmdE;
    // scattering applies directly to momentum (1st order)
    perpmomvar = scatvar*mom*mom;
  }

  template <class KTRAJ> double ElementXing<KTRAJ>::radiationFraction() const {
    double retval(0.0);
    for(auto const& mxing : matXings())
      retval += mxing.dmat_.radiationFraction(mxing.plen_);
    return retval;
  }

}
#endif
