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
      using STRAJPTR = std::shared_ptr<SensorLine>;
      // copy constructor
      ScintHit(ScintHit<KTRAJ> const& rhs):
          ResidualHit<KTRAJ>(rhs),
          tvar_(rhs.timeVariance()),
          wvar_(rhs.widthVariance()),
          tpca_(rhs.tpca_),
          rresid_(rhs.refResidual()){
        /**/
      };
      // clone op for reinstantiation
      std::shared_ptr< Hit<KTRAJ> > clone(CloneContext& context) const override{
        auto rv = std::make_shared< ScintHit<KTRAJ> >(*this);
        auto ca = rv->closestApproach();
        auto trajectory = std::make_shared<KTRAJ>(ca.particleTraj());
        ca.setTrajectory(trajectory);
        rv->setClosestApproach(ca);
        return rv;
      };
      // ResidualHit interface implementation
      unsigned nResid() const override { return 1; } // 1 time residual
      Residual const& refResidual(unsigned ires=0) const override;
      double time() const override { return tpca_.particleToca(); }
      void updateReference(PTRAJ const& ptraj) override;
      KTRAJPTR refTrajPtr() const override { return tpca_.particleTrajPtr(); }
      void updateState(MetaIterConfig const& config,bool first) override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
     // scintHit explicit interface
      ScintHit(PCA const& pca, double tvar, double wvar);
      virtual ~ScintHit(){}
      // the line encapsulates both the measurement value (through t0), and the light propagation model (through the velocity)
      auto sensorAxis() const { return tpca_.sensorTrajPtr(); }
      auto const& closestApproach() const { return tpca_; }
      double timeVariance() const { return tvar_; }
      double widthVariance() const { return wvar_; }
      auto precision() const { return tpca_.precision(); }
    private:
      double tvar_; // variance in the time measurement: assumed independent of propagation distance/time
      double wvar_; // variance in transverse position of the sensor/measurement in mm.  Assumes cylindrical error, could be more general
      CA tpca_; // reference time and position of closest approach to the axis
      Residual rresid_; // residual WRT most recent reference parameters

      // modifiers to support cloning
      void setClosestApproach(const CA& ca){ tpca_ = ca; }
  };

  template <class KTRAJ> ScintHit<KTRAJ>::ScintHit(PCA const& pca, double tvar, double wvar) :
    tvar_(tvar), wvar_(wvar),
    tpca_(static_cast<CA const&>(pca))
  {}

  template <class KTRAJ> Residual const& ScintHit<KTRAJ>::refResidual(unsigned ires) const {
    if(ires !=0)throw std::invalid_argument("Invalid residual");
    return rresid_;
  }

  template <class KTRAJ> void ScintHit<KTRAJ>::updateReference(PTRAJ const& ptraj) {
    // use previous hint, or initialize from the sensor time
    CAHint tphint = tpca_.usable() ?  tpca_.hint() : CAHint(sensorAxis()->measurementTime(), sensorAxis()->measurementTime());
    PCA pca(ptraj,sensorAxis(),tphint,precision());
    tpca_ = static_cast<CA>(pca);
    if(!tpca_.usable())throw std::runtime_error("ScintHit TPOCA failure");
  }

  template <class KTRAJ> void ScintHit<KTRAJ>::updateState(MetaIterConfig const& config,bool first) {
    // check that TPCA position is consistent with the physical sensor. This can be off if the CA algorithm finds the wrong helix branch
    // early in the fit when t0 has very large errors.
    // If it is unphysical try to adjust it back using a better hint.
    auto ppos = tpca_.particlePoca().Vect();
    auto const& sstart = sensorAxis()->start();
    auto const& send = sensorAxis()->end();
    // tolerance should come from the config.  Should also test relative to the error. FIXME
    double tol = sensorAxis()->length()*1.0;
    double sdist = (ppos - sensorAxis()->middle()).Dot(sensorAxis()->direction());
    if( (ppos-sstart).Dot(sensorAxis()->direction()) < -tol || (ppos-send).Dot(sensorAxis()->direction()) > tol) {
      // adjust hint to the middle and try agian
      double sspeed = tpca_.particleTraj().velocity(tpca_.particleToca()).Dot(sensorAxis()->direction());
      auto tphint = tpca_.hint();
      tphint.particleToca_ -= sdist/sspeed;
      tpca_ = CA(tpca_.particleTrajPtr(),sensorAxis(),tphint,precision());
      // should check if this is still unphysical and disable the hit if so FIXME
      sdist = (tpca_.particlePoca().Vect() - sensorAxis()->middle()).Dot(sensorAxis()->direction());
    }

    // residual is just delta-T at CA.
    // the variance includes the measurement variance and the tranvserse size (which couples to the relative direction)
    // Might want to do more updating (set activity) based on DOCA in future: TODO
    double dd2 = tpca_.dirDot()*tpca_.dirDot();
    double totvar = tvar_ + wvar_*dd2/(sensorAxis()->speed(sdist)*sensorAxis()->speed(sdist)*(1.0-dd2));
    rresid_ = Residual(tpca_.deltaT(),totvar,0.0,tpca_.dTdP());
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
      sensorAxis()->print(ost,detail);
    }
  }

}
#endif
