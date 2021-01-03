#ifndef KinKal_ScintHit_hh
#define KinKal_ScintHit_hh
//
//  simple example hit subclass representing a time measurement using scintillator light from a crystal or plastic scintillator
//
#include "KinKal/Detector/Hit.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/Detector/ElementXing.hh"
#include <stdexcept>
namespace KinKal {

  template <class KTRAJ> class ScintHit : public Hit<KTRAJ> {
    public:
      using HIT = Hit<KTRAJ>;
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using PTCA = PiecewiseClosestApproach<KTRAJ,Line>;
      using EXING = ElementXing<KTRAJ>;
      using EXINGPTR = std::shared_ptr<EXING>;

      // Hit interface overrrides
      unsigned nResid() const override { return 1; } // 1 time residual
      bool active(unsigned ires) const override;
      Residual residual(unsigned ires=0) const override;
      double time() const override { return tpoca_.particleToca(); }
      void update(PKTRAJ const& pktraj) override;
      void update(PKTRAJ const& pktraj, MetaIterConfig const& config) override;
      EXINGPTR const& detXingPtr() const override { return null_; }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // scintHit explicit interface
      ScintHit(Line const& sensorAxis, double tvar, double wvar) : 
	saxis_(sensorAxis), tvar_(tvar), wvar_(wvar), active_(true), precision_(1e-6) {}
      virtual ~ScintHit(){}
      Residual const& timeResidual() const { return rresid_;
    // the line encapsulates both the measurement value (through t0), and the light propagation model (through the velocity)
      Line const& sensorAxis() const { return saxis_; }
      PTCA const& closestApproach() const { return tpoca_; }
      double timeVariance() const { return tvar_; }
      double widthVariance() const { return wvar_; }
    private:
      Line saxis_; // symmetry axis of this sensor
      double tvar_; // variance in the time measurement: assumed independent of propagation distance/time
      double wvar_; // variance in transverse position of the sensor/measurement in mm.  Assumes cylindrical error, could be more general
      bool active_; // active or not
      EXINGPTR null_; // no detector material xing: should be added
      PTCA tpoca_; // reference time and distance of closest approach to the axis.
      // caches
      Residual rresid_; // residual WRT most recent reference parameters
      double precision_; // current precision
  };

  template <class KTRAJ> bool ScintHit<KTRAJ>::active(unsigned ires) const {
    if(ires == 0 && active_)
      return true;
    else
      return false;
  }

  template <class KTRAJ> Residual const& ScintHit<KTRAJ>::residual(unsigned ires) const {
    if(ires !=0)throw std::invalid_argument("Invalid residual");
    return rresid_;
  }

  template <class KTRAJ> void ScintHit<KTRAJ>::update(PKTRAJ const& pktraj) {
    // compute PTCA
    CAHint tphint( saxis_.t0(), saxis_.t0());
    if(tpoca_.usable()) tphint = CAHint(tpoca_.particleToca(),tpoca_.sensorToca());
    tpoca_ = PTCA(pktraj,saxis_,tphint,precision_);
    if(tpoca_.usable()){
      // residual is just delta-T at CA. 
      // the variance includes the measurement variance and the tranvserse size (which couples to the relative direction)
      double dd2 = tpoca.dirDot()*tpoca.dirDot();
      double totvar = tvar_ + wvar_*dd2/(saxis_.speed()*saxis_.speed()*(1.0-dd2));
      rresid_ = Residual(Residual::dtime,tpoca.tpData(),tpoca.deltaT(),totvar,-tpoca.dTdP());
      reftraj_ = pktraj.nearestPiece(tpoca_.particleToca());
    } else
      throw std::runtime_error("PTCA failure");
  }

  template <class KTRAJ> void ScintHit<KTRAJ>::update(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) {
    // for now, no updates are needed.  Eventually could test for consistency, update errors, etc
    precision_ = miconfig.tprec_;
    update(pktraj);
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
