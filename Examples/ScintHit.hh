#ifndef KinKal_ScintHit_hh
#define KinKal_ScintHit_hh
//
//  simple hit subclass representing a time measurement using scintillator light from a crystal or plastic scintillator
//
#include "KinKal/Detector/ResidualHit.hh"
#include "KinKal/Trajectory/SensorLine.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include <stdexcept>
namespace KinKal {

  template <class KTRAJ> class ScintHit : public ResidualHit<KTRAJ> {
    public:
      using PTRAJ = ParticleTrajectory<KTRAJ>;
      using PCA = PiecewiseClosestApproach<KTRAJ,SensorLine>;
      using CA = ClosestApproach<KTRAJ,SensorLine>;
      using RESIDHIT = ResidualHit<KTRAJ>;
      using HIT = Hit<KTRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      // ResidualHit interface implementation
      unsigned nResid() const override { return 1; } // 1 time residual
      Residual const& refResidual(unsigned ires=0) const override;
      double time() const override { return tpca_.particleToca(); }
      void updateReference(KTRAJPTR const& ktrajptr) override;
      KTRAJPTR const& refTrajPtr() const override { return tpca_.particleTrajPtr(); }
      void updateState(MetaIterConfig const& config,bool first) override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
     // scintHit explicit interface
      ScintHit(PCA const& pca, double tvar, double wvar);
      virtual ~ScintHit(){}
      // the line encapsulates both the measurement value (through t0), and the light propagation model (through the velocity)
      auto const& sensorAxis() const { return saxis_; }
      auto const& closestApproach() const { return tpca_; }
      double timeVariance() const { return tvar_; }
      double widthVariance() const { return wvar_; }
      auto precision() const { return tpca_.precision(); }
    private:
      SensorLine saxis_; // symmetry axis of this sensor
      double tvar_; // variance in the time measurement: assumed independent of propagation distance/time
      double wvar_; // variance in transverse position of the sensor/measurement in mm.  Assumes cylindrical error, could be more general
      CA tpca_; // reference time and position of closest approach to the axis
      Residual rresid_; // residual WRT most recent reference parameters
  };

  template <class KTRAJ> ScintHit<KTRAJ>::ScintHit(PCA const& pca, double tvar, double wvar) :
    saxis_(pca.sensorTraj()), tvar_(tvar), wvar_(wvar),
    tpca_(pca.localTraj(),saxis_,pca.precision(),pca.tpData(),pca.dDdP(),pca.dTdP())
  {}

  template <class KTRAJ> Residual const& ScintHit<KTRAJ>::refResidual(unsigned ires) const {
    if(ires !=0)throw std::invalid_argument("Invalid residual");
    return rresid_;
  }

  template <class KTRAJ> void ScintHit<KTRAJ>::updateReference(KTRAJPTR const& ktrajptr) {
    // use previous hint, or initialize from the sensor time
    CAHint tphint = tpca_.usable() ?  tpca_.hint() : CAHint(saxis_.measurementTime(), saxis_.measurementTime());
    tpca_ = CA(ktrajptr,saxis_,tphint,precision());
    if(!tpca_.usable())throw std::runtime_error("ScintHit TPOCA failure");
  }

  template <class KTRAJ> void ScintHit<KTRAJ>::updateState(MetaIterConfig const& config,bool first) {
    // check that TPCA position is consistent with the physical sensor. This can be off if the CA algorithm finds the wrong helix branch
    // early in the fit when t0 has very large errors.
    // If it is unphysical try to adjust it back using a better hint.
    auto ppos = tpca_.particlePoca().Vect();
    auto const& sstart = saxis_.start();
    auto const& send = saxis_.end();
    // tolerance should come from the config.  Should also test relative to the error. FIXME
    double tol = saxis_.length()*1.0;
    double sdist = (ppos - saxis_.middle()).Dot(saxis_.direction());
    if( (ppos-sstart).Dot(saxis_.direction()) < -tol || (ppos-send).Dot(saxis_.direction()) > tol) {
      // adjust hint to the middle and try agian
      double sspeed = tpca_.particleTraj().velocity(tpca_.particleToca()).Dot(saxis_.direction());
      auto tphint = tpca_.hint();
      tphint.particleToca_ -= sdist/sspeed;
      tpca_ = CA(tpca_.particleTrajPtr(),saxis_,tphint,precision());
      // should check if this is still unphysical and disable the hit if so FIXME
      sdist = (tpca_.particlePoca().Vect() - saxis_.middle()).Dot(saxis_.direction());
    }

    // residual is just delta-T at CA.
    // the variance includes the measurement variance and the tranvserse size (which couples to the relative direction)
    // Might want to do more updating (set activity) based on DOCA in future: TODO
    double dd2 = tpca_.dirDot()*tpca_.dirDot();
    double totvar = tvar_ + wvar_*dd2/(saxis_.speed(sdist)*saxis_.speed(sdist)*(1.0-dd2));
    rresid_ = Residual(tpca_.deltaT(),totvar,0.0,true,tpca_.dTdP());
    this->updateWeight(config);
  }

  template<class KTRAJ> void ScintHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    if(this->active())
      ost<<"Active ";
    else
      ost<<"Inactive ";
    ost << " ScintHit  tvar " << tvar_ << " wvar " << wvar_ << std::endl;
    if(detail > 0){
      ost << "Line ";
      saxis_.print(ost,detail);
    }
  }

}
#endif
