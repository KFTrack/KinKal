#ifndef KinKal_ParameterHit_hh
#define KinKal_ParameterHit_hh
//
//  direct constraint on a subset of parameters expressed as a 'hit'.  This allows
//  external information to be added to the fit.
//
#include "KinKal/Detector/Hit.hh"
#include "KinKal/General/Vectors.hh"
//#include "KinKal/General/Parameters.hh"
#include <stdexcept>
namespace KinKal {

  template <class KTRAJ> class ParameterHit : public Hit<KTRAJ> {
    public:
      using PMASK = std::array<bool,NParams()>; // parameter mask
      using HIT = Hit<KTRAJ>;
      using PTRAJ = ParticleTrajectory<KTRAJ>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;

      // Hit interface overrrides
      bool active() const override { return ncons_ > 0; }
      Chisq chisq(Parameters const& pdata) const override;
      double time() const override { return time_; }
      void updateState(MetaIterConfig const& config,bool first) override;
      Weights const& weight() const override { return weight_; }
      // parameter constraints are absolute and can't be updated
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      void updateReference(KTRAJPTR const& ktrajptr) override { reftraj_ = ktrajptr; }
      KTRAJPTR const& refTrajPtr() const override { return reftraj_; }
     // ParameterHit-specfic interface
      // construct from constraint values, time, and mask of which parameters to constrain
      ParameterHit(double time, PTRAJ const& ptraj, Parameters const& params, PMASK const& pmask);
      virtual ~ParameterHit(){}
      unsigned nDOF() const override { return ncons_; }
      Parameters const& constraintParameters() const { return params_; }
      PMASK const& constraintMask() const { return pmask_; }
    private:
      double time_; // time of this constraint: must be supplied on construction and does not change
      KTRAJPTR reftraj_; // reference WRT this hits weight was calculated
      Parameters params_; // constraint parameters with covariance
      Weights pweight_; // weight from these (masked) parameters
      Weights weight_; // current weight, including temp effects
      PMASK pmask_; // subset of parmeters to constrain
      DMAT mask_; // matrix to mask off unconstrainted parameters
      unsigned ncons_; // number of parameters constrained
  };

  template<class KTRAJ> ParameterHit<KTRAJ>::ParameterHit(double time, PTRAJ const& ptraj, Parameters const& params, PMASK const& pmask) :
    time_(time), reftraj_(ptraj.nearestTraj(time)), params_(params), pmask_(pmask), ncons_(0) {
      // create the mask matrix; Use a temporary, not the data member, as root has caching problems with that (??)
      mask_ = ROOT::Math::SMatrixIdentity();
      // count constrained parameters, and mask off unused parameters
      for(size_t ipar=0;ipar < NParams(); ipar++){
        if(pmask_[ipar]){
          ncons_++;
        } else {
          mask_(ipar,ipar) = 0.0;
        }
      }
      // Mask Off unused parameters
      // 2 steps needed here, as otherwise root caching results in incomplete objects
      Weights weight(params);
      DMAT wmat = weight.weightMat();
      wmat = ROOT::Math::Similarity(mask_,wmat);
      DVEC wvec = weight.weightVec();
      DVEC wreduced = wvec*mask_;
      pweight_ = Weights(wreduced, wmat);
    }

  template <class KTRAJ> void ParameterHit<KTRAJ>::updateState(MetaIterConfig const& miconfig,bool first) {
    weight_ = pweight_; // do this in 2 steps to avoid SMatrix caching issue
    weight_ *= 1.0/miconfig.varianceScale();
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
    if(!wmat.Invert()) throw std::runtime_error("ParameterHit inversion failure");
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
