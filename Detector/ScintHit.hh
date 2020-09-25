
#ifndef KinKal_ScintHit_hh
#define KinKal_ScintHit_hh
//
//  class representing a timing measurement using scintillator light from a crystal or plastic scintillator
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/Detector/DetectorHit.hh"
#include "KinKal/Trajectory/Line.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "KinKal/Detector/DetectorXing.hh"
#include <stdexcept>
namespace KinKal {

  template <class KTRAJ> class ScintHit : public DetectorHit<KTRAJ> {
    public:
      using DHIT = DetectorHit<KTRAJ>;
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using PTCA = PiecewiseClosestApproach<KTRAJ,Line>;
      using DXING = DetectorXing<KTRAJ>;
      using DXINGPTR = std::shared_ptr<DXING>;

      // DetectorHit interface overrrides
      Weights weight() const override;
      double chi(Parameters const& pdata) const override;
      double time() const override { return rresid_.time(); }
      void update(PKTRAJ const& pktraj, MetaIterConfig const& config) override;
      void update(PKTRAJ const& pktraj) override;
      bool isActive() const override { return active_; }
      DXINGPTR const& detXingPtr() const override { return null_; }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      unsigned nDOF() const override { return 1; }
      // the line encapsulates both the measurement value (through t0), and the light propagation model (through the velocity)
      Line const& sensorAxis() const { return saxis_; }
      ScintHit(Line const& sensorAxis, double tvar, double wvar) : 
	saxis_(sensorAxis), tvar_(tvar), wvar_(wvar), active_(true), precision_(1e-6) {}
      virtual ~ScintHit(){}
      double timeVariance() const { return tvar_; }
      double widthVariance() const { return wvar_; }
      Residual const& refResidual() const { return rresid_; }
      Parameters const& refParams() const { return rparams_; }
    private:
      Line saxis_;
      double tvar_; // variance in the time measurement: assumed independent of propagation distance/time FIXME!
      double wvar_; // variance in transverse position of the sensor/measurement in mm.  Assumes cylindrical error, Should be more general FIXME!
      bool active_; // active or not (pat. rec. tool)
      DXINGPTR null_; // no detector material xing: should be added TODO
      // caches
      Residual rresid_; // residual WRT most recent reference parameters
      Parameters rparams_; // reference parameters
      double precision_; // current precision
  };

  template <class KTRAJ> Weights ScintHit<KTRAJ>::weight() const {
    if(isActive()){
      // convert derivatives to a Nx1 matrix (for root)
      ROOT::Math::SMatrix<double,NParams(),1> dRdPM;
      dRdPM.Place_in_col(rresid_.dRdP(),0,0);
      // convert the variance into a 1X1 matrix
      ROOT::Math::SMatrix<double, 1,1, ROOT::Math::MatRepSym<double,1>> RVarM;
      // weight by inverse variance
      double tvar = rresid_.variance();
      RVarM(0,0) = 1.0/tvar;
      // expand these into the weight matrix
      DMAT wmat = ROOT::Math::Similarity(dRdPM,RVarM);
      // translate residual value into weight vector WRT the reference parameters
      // sign convention reflects resid = measurement - prediction
      DVEC wvec = wmat*rparams_.parameters() + rresid_.dRdP()*rresid_.value()/tvar;
      return Weights(wvec,wmat);
   } else
      return Weights();
  }

  template <class KTRAJ> void ScintHit<KTRAJ>::update(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) {
    // for now, no updates are needed.  Eventually could be tests for consistency TODO
    precision_ = miconfig.tprec_;
    update(pktraj);
  }

  template <class KTRAJ> void ScintHit<KTRAJ>::update(PKTRAJ const& pktraj) {
      // compute PTCA
      CAHint tphint( saxis_.t0(), saxis_.t0());
      PTCA tpoca(pktraj,saxis_,tphint,precision_);
      if(tpoca.usable()){
	// residual is just delta-T at CA. 
	// the variance includes the measurement variance and the tranvserse size (which couples to the relative direction)
	double dd2 = tpoca.dirDot()*tpoca.dirDot();
	double totvar = tvar_ + wvar_*dd2/(saxis_.speed()*saxis_.speed()*(1.0-dd2));
	rresid_ = Residual(Residual::dtime,tpoca.tpData(),tpoca.deltaT(),totvar,-tpoca.dTdP());
	rparams_ = pktraj.nearestPiece(rresid_.time()).params();
      } else
	throw std::runtime_error("CA failure");
  }

  template <class KTRAJ> double ScintHit<KTRAJ>::chi(Parameters const& pdata) const {
    double retval(0.0);
    if(isActive()) {
      // compute the difference between these parameters and the reference parameters
      DVEC dpvec = pdata.parameters() - rparams_.parameters(); 
      // use the differnce to 'correct' the reference residual to be WRT these parameters
      double uresid = rresid_.value() - ROOT::Math::Dot(dpvec,rresid_.dRdP());
      // project the parameter covariance into a residual space variance
      double rvar = ROOT::Math::Similarity(rresid_.dRdP(),pdata.covariance());
      // add the measurement variance
      rvar +=  rresid_.variance();
      // chi is the ratio of these
      retval = uresid/sqrt(rvar);
    }
    return retval;
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
