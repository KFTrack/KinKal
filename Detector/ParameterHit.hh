#ifndef KinKal_ParameterHit_hh
#define KinKal_ParameterHit_hh
//
//  direct constraint on a subset of parameters expressed as a 'hit'.  This allows
//  external information to be added to the fit.
//
#include "KinKal/Detector/Hit.hh"
#include "KinKal/General/Vectors.hh"
#include <stdexcept>
namespace KinKal {

  template <class KTRAJ> class ParameterHit : public Hit<KTRAJ> {
    public:
      using PMASK = std::array<bool,NParams()>; // parameter mask
      using HIT = Hit<KTRAJ>;
      using PKTRAJ = ParticleTrajectory<KTRAJ>;

      // Hit interface overrrides
      bool active() const override { return ncons_ > 0; }
      Chisq chisq(Parameters const& pdata) const override;
      double time() const override { return time_; }
      void update(PKTRAJ const& pktraj) override { HIT::refparams_ = pktraj.nearestPiece(time()).params(); }
      void update(PKTRAJ const& pktraj, MetaIterConfig const& miconfig) override {
        HIT::wscale_ = 1.0/miconfig.varianceScale(); update(pktraj); }
      // parameter constraints are absolute and can't be updated
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // ParameterHit-specfic interface
      // construct from constraint values, time, and mask of which parameters to constrain
      ParameterHit(double time, Parameters const& params, PMASK const& pmask);
      virtual ~ParameterHit(){}
      unsigned nDOF() const { return ncons_; }
      Parameters const& constraintParameters() const { return params_; }
      PMASK const& constraintMask() const { return pmask_; }
    private:
      double time_; // time of this constraint: must be supplied on construction and does not change
      Parameters params_; // constraint parameters with covariance
      PMASK pmask_; // subset of parmeters to constrain
      DMAT mask_; // matrix to mask off unconstrainted parameters
      unsigned ncons_; // number of parameters constrained
  };

  template<class KTRAJ> ParameterHit<KTRAJ>::ParameterHit(double time, Parameters const& params, PMASK const& pmask) :
    time_(time), params_(params), pmask_(pmask), mask_(ROOT::Math::SMatrixIdentity()), ncons_(0) {

      Weights weight(params);
      // count constrained parameters, and mask off unused parameters
      for(size_t ipar=0;ipar < NParams(); ipar++){
        if(pmask_[ipar]){
          ncons_++;
        } else {
          // zero unconstrained values
          mask_(ipar,ipar) = 0.0;
        }
      }
      // Mask Off unused parameters
      DMAT wmat = ROOT::Math::Similarity(mask_,weight.weightMat());
      // 2 steps needed her, as otherwise root caching results in incomplete objects
      DVEC wvec = weight.weightVec();
      DVEC wreduced = wvec*mask_;
      HIT::weight_ = Weights(wreduced, wmat);
    }

  template <class KTRAJ> Chisq ParameterHit<KTRAJ>::chisq(Parameters const& pdata) const {
    // chi measures the dimensionless tension between this constraint and the given parameters, including uncertainty
    // on both the measurement and the trajectory estimate.
    // Compute as the parameter difference contracted through the sum covariance of the 2.
    Parameters pdiff = pdata;
    pdiff.parameters() *= -1.0; // so I can subtract in the next step
    pdiff += params_;  // this is now the difference of parameters but sum of covariances
    // invert the covariance matrix
    DMAT wmat = pdiff.covariance();
    if(!wmat.Invert()) throw std::runtime_error("Inversion failure");
    // zero out unconstrainted parts
    wmat = ROOT::Math::Similarity(mask_,wmat);
    double chisq = ROOT::Math::Similarity(pdiff.parameters(),wmat);
    // sign of chi has no meaning here
    return Chisq(chisq,nDOF());
  }

  template<class KTRAJ> void ParameterHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    if(this->active())
      ost<<"Active ";
    else
      ost<<"Inactive ";
    ost << " ParameterHit Hit" << std::endl;
    if(detail > 0){
      for(size_t ipar=0;ipar < NParams(); ipar++){
        auto tpar = static_cast<typename KTRAJ::ParamIndex>(ipar);
        if (pmask_[ipar]) {
          ost << " constraining parameter " << KTRAJ::paramName(tpar) << " to value " << params_.parameters()[ipar] << " +- " << sqrt(params_.covariance()(ipar,ipar)) << std::endl;
        }
      }
    }
  }

} // namespace KinKal
#endif
