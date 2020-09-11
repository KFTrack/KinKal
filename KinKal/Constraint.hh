#ifndef KinKal_Constraint_hh
#define KinKal_Constraint_hh
//
//  direct constraint on a subset of parameters expressed as a 'measurement'
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/DetectorHit.hh"
#include "KinKal/Vectors.hh"
#include <stdexcept>
namespace KinKal {

  template <class KTRAJ> class Constraint : public DetectorHit<KTRAJ> {
    public:
      using PMASK = std::array<bool,NParams()>; // parameter mask
      using DHIT = DetectorHit<KTRAJ>;
      using PKTRAJ = ParticleTrajectory<KTRAJ>;
      using DXING = DetectorXing<KTRAJ>;
      using DXINGPTR = std::shared_ptr<DXING>;

      // Measurement interface overrrides
      Weights weight() const override { return weight_; }
      unsigned nDOF() const override { return ncons_; }
      double chi(Parameters const& pdata) const override;
      double time() const override { return time_; }
      // constraints are absolute and can't be updatedA
      void update(PKTRAJ const& pktraj, MetaIterConfig const& config) override {}
      void update(PKTRAJ const& pktraj) override {}
      bool isActive() const override { return ncons_ > 0; }
      DXINGPTR const& detXingPtr() const override { return null_; }
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      // construct from constraint values, time, and mask of which parameters to constrain
      Constraint(double time, Parameters const& params, PMASK const& pmask);
      virtual ~Constraint(){}
      // Constraint-specfic interface
      Parameters const& constraintParameters() const { return params_; }
      PMASK const& constraintMask() const { return pmask_; }
    private:
      double time_; // time of this constraint: must be supplied on construction
      Parameters params_; // constraint parameters with covariance
      Weights weight_; // constraint WRT reference parameters expressed as a weight 
      PMASK pmask_; // subset of parmeters to constrain
      DMAT mask_; // matrix to mask of unconstrainted parameters
      unsigned ncons_; // number of parameters constrained
      DXINGPTR null_; // null detector material crossing 
  };

  template<class KTRAJ> Constraint<KTRAJ>::Constraint(double time, Parameters const& params, PMASK const& pmask) :
    time_(time), params_(params), pmask_(pmask), mask_(ROOT::Math::SMatrixIdentity()), ncons_(0) {
      weight_ = Weights(params_);
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
      weight_.weightMat() = ROOT::Math::Similarity(mask_,weight_.weightMat());
      weight_.weightVec() = weight_.weightVec()*mask_;
    }

  template <class KTRAJ> double Constraint<KTRAJ>::chi(Parameters const& pdata) const {
  // chi measures the dimensionless 'distance' between this constraint and the given parameters.
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
    return sqrt(chisq);
  }

  template<class KTRAJ> void Constraint<KTRAJ>::print(std::ostream& ost, int detail) const {
    if(this->isActive())
      ost<<"Active ";
    else
      ost<<"Inactive ";
    ost << " Constraint Hit" << std::endl;
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
