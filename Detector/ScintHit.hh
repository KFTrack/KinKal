#ifndef KinKal_ScintHit_hh
#define KinKal_ScintHit_hh
//
//  simple hit subclass representing a time measurement using scintillator light from a crystal or plastic scintillator
//
#include "KinKal/Detector/ResidualHit.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include <stdexcept>
namespace KinKal {

  template <class KTRAJ> class ScintHit : public ResidualHit<KTRAJ> {
    public:
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using PCA = PiecewiseClosestApproach<KTRAJ,Line>;
      using CA = ClosestApproach<KTRAJ,Line>;
      using RESIDHIT = ResidualHit<KTRAJ>;
      using HIT = Hit<KTRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      // Hit interface implementation
      unsigned nResid() const override { return 1; } // 1 time residual
      bool activeRes(unsigned ires=0) const override;
      Residual const& residual(unsigned ires=0) const override;
      double time() const override { return tpca_.particleToca(); }
      void updateReference(KTRAJPTR const& ktrajptr) override;
      KTRAJPTR const& refTrajPtr() const override { return tpca_.particleTrajPtr(); }
      void updateState(MetaIterConfig const& config,bool first) override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // scintHit explicit interface
      ScintHit(PCA const& pca, double tvar, double wvar);
      virtual ~ScintHit(){}
      auto const& timeResidual() const { return rresid_; }
      // the line encapsulates both the measurement value (through t0), and the light propagation model (through the velocity)
      auto const& sensorAxis() const { return saxis_; }
      auto const& closestApproach() const { return tpca_; }
      double timeVariance() const { return tvar_; }
      double widthVariance() const { return wvar_; }
      auto precision() const { return tpca_.precision(); }
    private:
      Line saxis_; // symmetry axis of this sensor
      double tvar_; // variance in the time measurement: assumed independent of propagation distance/time
      double wvar_; // variance in transverse position of the sensor/measurement in mm.  Assumes cylindrical error, could be more general
      bool active_; // active or not
      CA tpca_; // reference time and position of closest approach to the axis
      Residual rresid_; // residual WRT most recent reference parameters
  };

  template <class KTRAJ> ScintHit<KTRAJ>::ScintHit(PCA const& pca, double tvar, double wvar) :
    saxis_(pca.sensorTraj()), tvar_(tvar), wvar_(wvar), active_(true),
    tpca_(pca.localTraj(),saxis_,pca.precision(),pca.tpData(),pca.dDdP(),pca.dTdP())
  {}

  template <class KTRAJ> bool ScintHit<KTRAJ>::activeRes(unsigned ires) const {
    if(ires == 0 && active_)
      return true;
    else
      return false;
  }

  template <class KTRAJ> Residual const& ScintHit<KTRAJ>::residual(unsigned ires) const {
    if(ires !=0)throw std::invalid_argument("Invalid residual");
    return rresid_;
  }

  template <class KTRAJ> void ScintHit<KTRAJ>::updateReference(KTRAJPTR const& ktrajptr) {
    // compute PCA
    CAHint tphint( saxis_.t0(), saxis_.t0());
    // don't update the hint: initial T0 values can be very poor, which can push the CA calculation onto the wrong helix loop,
    // from which it's impossible to ever get back to the correct one.  Active loop checking might be useful eventually too TODO
    //    if(tpca_.usable()) tphint = CAHint(tpca_.particleToca(),tpca_.sensorToca());
    tpca_ = CA(ktrajptr,saxis_,tphint,precision());
    if(!tpca_.usable())throw std::runtime_error("ScintHit TPOCA failure");
    }

  template <class KTRAJ> void ScintHit<KTRAJ>::updateState(MetaIterConfig const& config,bool first) {
    // residual is just delta-T at CA.
    // the variance includes the measurement variance and the tranvserse size (which couples to the relative direction)
    // Might want to do more updating (set activity) based on DOCA in future: TODO
    double dd2 = tpca_.dirDot()*tpca_.dirDot();
    double totvar = tvar_ + wvar_*dd2/(saxis_.speed()*saxis_.speed()*(1.0-dd2));
    rresid_ = Residual(tpca_.deltaT(),totvar,-tpca_.dTdP());
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
