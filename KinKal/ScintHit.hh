
#ifndef KinKal_ScintHit_hh
#define KinKal_ScintHit_hh
//
//  class representing a timing measurement using scintillator light from a crystal or plastic scintillator
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/DetectorHit.hh"
#include "KinKal/DistanceToTime.hh"
#include "KinKal/Line.hh"
#include "KinKal/PieceClosestApproach.hh"
#include <stdexcept>
namespace KinKal {

  template <class KTRAJ> class ScintHit : public DetectorHit<KTRAJ> {
    public:
      using THIT = DetectorHit<KTRAJ>;
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using PTPOCA = PieceClosestApproach<KTRAJ,Line>;
      
      // Hit interface overrrides
      virtual void resid(PKTRAJ const& pktraj, Residual& resid,double precision) const override;
      virtual void update(PKTRAJ const& pktraj, MetaIterConfig const& config, Residual& resid) override;
      virtual unsigned nDOF() const override { return 1; }
//      virtual double tension() const override { return tpoca_.doca()/sqrt(wvar_); } 
      virtual double tension() const override { return 0.0; }  // FIXME!
      virtual void print(std::ostream& ost=std::cout,int detail=0) const override;
      // the line encapsulates both the measurement value (through t0), and the light propagation model (through the velocity)
      Line const& sensorAxis() const { return saxis_; }
      ScintHit(Line const& sensorAxis, double tvar, double wvar, bool active=true) : 
	THIT(active), saxis_(sensorAxis), tvar_(tvar), wvar_(wvar) {}
      virtual ~ScintHit(){}
      double timeVariance() const { return tvar_; }
      double widthVariance() const { return wvar_; }
    private:
      Line saxis_;
      double tvar_; // variance in the time measurement: assumed independent of propagation distance/time FIXME!
      double wvar_; // variance in transverse position of the sensor/measurement in mm.  Assumes cylindrical error, Should be more general FIXME!
  };

  template <class KTRAJ> void ScintHit<KTRAJ>::resid(PKTRAJ const& pktraj,  Residual& resid, double precision) const {
    // compute PTPOCA
    CAHint tphint( saxis_.t0(), saxis_.t0());
    PTPOCA tpoca(pktraj,saxis_,tphint,precision);
 
    if(tpoca.usable()){
      // residual is just delta-T at POCA. 
      // the variance includes the measurement variance and the tranvserse size (which couples to the relative direction)
	double dd2 = tpoca.dirDot()*tpoca.dirDot();
	double totvar = tvar_ + wvar_*dd2/(saxis_.speed()*saxis_.speed()*(1.0-dd2));
	resid = Residual(Residual::dtime,tpoca.tpData(),tpoca.deltaT(),totvar,-tpoca.dTdP());
    } else
      throw std::runtime_error("POCA failure");
  }

  template <class KTRAJ> void ScintHit<KTRAJ>::update(PKTRAJ const& pktraj, MetaIterConfig const& miconfig, Residual& residual) {
  // for now, no updates are needed.  Eventually could be tests for consistency, tension etc FIXME!
    resid(pktraj,residual,miconfig.tprec_);
  }

  template<class KTRAJ> void ScintHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    if(this->isActive())
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
