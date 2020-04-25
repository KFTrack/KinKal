
#ifndef KinKal_ScintHit_hh
#define KinKal_ScintHit_hh
//
//  class representing a timing measurement using scintillator light from a crystal or plastic scintillator
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/THit.hh"
#include "KinKal/D2T.hh"
#include "KinKal/TLine.hh"
#include "KinKal/TPoca.hh"
#include "KinKal/LRAmbig.hh"
#include <stdexcept>
namespace KinKal {

  template <class KTRAJ> class ScintHit : public THit<KTRAJ> {
    public:
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef Residual<KTRAJ::NParams()> RESIDUAL;
      typedef TPoca<PKTRAJ,TLine> TPOCA;
      typedef typename KTRAJ::DVEC DVEC; 
      // THit interface overrrides
      virtual void resid(PKTRAJ const& pktraj, RESIDUAL& resid) const override;
      virtual void update(PKTRAJ const& pktraj, MConfig const& config, RESIDUAL& resid) override;
      virtual unsigned nDOF() const override { return 1; }
//      virtual float tension() const override { return tpoca_.doca()/sqrt(wvar_); } 
      virtual float tension() const override { return 0.0; }  // FIXME!
      virtual void print(std::ostream& ost=std::cout,int detail=0) const override;
      // the line encapsulates both the measurement value (through t0), and the light propagation model (through the velocity)
      TLine const& sensorAxis() const { return saxis_; }
      ScintHit(TLine const& sensorAxis, float tvar, float wvar, bool active=true) : 
	THit<KTRAJ>(active), saxis_(sensorAxis), tvar_(tvar), wvar_(wvar) {}
      virtual ~ScintHit(){}
      float timeVariance() const { return tvar_; }
      float widthVariance() const { return wvar_; }
    private:
      TLine saxis_;
      float tvar_; // variance in the time measurement: assumed independent of propagation distance/time FIXME!
      float wvar_; // variance in transverse position of the sensor/measurement in mm.  Assumes cylindrical error, Should be more general FIXME!
  };

  template <class KTRAJ> void ScintHit<KTRAJ>::resid(PKTRAJ const& pktraj,  RESIDUAL& resid) const {
    // compute TPOCA
    TPocaHint tphint;
    tphint.particleHint_ = true;
    tphint.particleToca_ = saxis_.t0();
    TPOCA tpoca(pktraj,saxis_,tphint);
 
    if(tpoca.usable()){
      // residual is just delta-T at POCA. 
      // the variance includes the measurement variance and the tranvserse size (which couples to the relative direction)
	double dd2 = tpoca.dirDot()*tpoca.dirDot();
	double totvar = tvar_ + wvar_*dd2/(saxis_.speed()*saxis_.speed()*(1.0-dd2));
	resid = RESIDUAL(RESIDUAL::dtime,tpoca,tpoca.deltaT(),totvar,tpoca.dTdP());
    } else
      throw std::runtime_error("POCA failure");
  }

  template <class KTRAJ> void ScintHit<KTRAJ>::update(PKTRAJ const& pktraj, MConfig const& mconfig, RESIDUAL& residual) {
  // for now, no updates are needed.  Eventually could be tests for consistency, tension etc FIXME!
    resid(pktraj,residual);
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
