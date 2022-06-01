#ifndef KinKal_ResidualHit_hh
#define KinKal_ResidualHit_hh
//
//  class representing a hit based on a set of uncorrelated but related (same sensor) measurements, expressed as residuals.
//  Each residual esimates the 1-dimenaional tension between the measurement and the value predicted by the reference trajectory
//
#include "KinKal/Detector/Hit.hh"
#include "KinKal/Detector/Residual.hh"
namespace KinKal {

  template <class KTRAJ> class ResidualHit : public Hit<KTRAJ> {
    public:
      // override of some Hit interface.  Subclasses must still implement update and material methods
      using HIT = Hit<KTRAJ>;
      using PTRAJ = ParticleTrajectory<KTRAJ>;
      bool active() const override { return nDOF() > 0; }
      Chisq chisq(Parameters const& params) const override;
      void updateWeight(MetaIterConfig const& config) override;
      Weights const& weight() const override { return weight_; }
      // ResidualHit specific interface.
      unsigned nDOF() const;
      // describe residuals associated with this hit
      virtual unsigned nResid() const = 0;
      // individual residuals may be active or inactive
      virtual bool activeRes(unsigned ires) const = 0;
      // reference residuals for this hit.  iDOF indexs the measurement and is hit-specific, outside the range will throw
      // this is generally biased as the refefence includes the effect of this hit
      virtual Residual const& residual(unsigned ires) const = 0;
      // residuals corrected to refer to the given set of parameters (1st-order)
      Residual residual(Parameters const& params, unsigned ires) const;
      // unbiased residuals WRT the reference parameters
      Residual unbiasedResidual(unsigned ires) const;
      // unbiased pull of this residual (including the uncertainty on the reference parameters)
      double pull(unsigned ires) const;
      ResidualHit() {}
    private:
      Weights weight_; // weight of this hit computed from the residuals
  };

  template <class KTRAJ> Residual ResidualHit<KTRAJ>::residual(Parameters const& pdata,unsigned ires) const {
    auto const& resid = residual(ires);
    // compute the difference between these parameters and the reference parameters
    DVEC dpvec = pdata.parameters() - HIT::referenceParameters().parameters();
    // project the parameter differnce to residual space and 'correct' the reference residual to be WRT these parameters
    double uresid = resid.value() - ROOT::Math::Dot(dpvec,resid.dRdP());
    return Residual(uresid,resid.variance(),resid.dRdP());
  }

  template <class KTRAJ> Residual ResidualHit<KTRAJ>::unbiasedResidual(unsigned ires) const {
    return residual(HIT::unbiasedParameters(),ires);
  }

  template <class KTRAJ> double ResidualHit<KTRAJ>::pull(unsigned ires) const {
    auto uparams = HIT::unbiasedParameters();
    auto ures = residual(uparams,ires);
    // project the parameter covariance into residual space
    double pvar = ROOT::Math::Similarity(ures.dRdP(),uparams.covariance());
    double pval = ures.value()/sqrt(ures.variance() + pvar);
    return pval;
  }

  template <class KTRAJ> Chisq ResidualHit<KTRAJ>::chisq(Parameters const& params) const {
    double chisq(0.0);
    unsigned ndof(0);
    for(unsigned ires=0; ires< nResid(); ires++) {
      if(activeRes(ires)) {
        ndof++;
        // compute the residual WRT the given parameters
        auto res = residual(params,ires);
        // project the parameter covariance into a residual space variance
        double rvar = ROOT::Math::Similarity(res.dRdP(),params.covariance());
        // check for unphysical values
        if(rvar<0){
//          std::cout << "neg resid var " << rvar << std::endl;
//          rvar = 0.0;
//          chisq = -1.0;
//          break;
          throw std::runtime_error("Covariance projection inconsistency");
        }
        // add the measurement variance
        rvar +=  res.variance();
        // add chisq for this DOF
        chisq += (res.value()*res.value())/rvar;
      }
    }
    return Chisq(chisq,ndof);
  }

  template <class KTRAJ> unsigned ResidualHit<KTRAJ>::nDOF() const {
    // each residual is counted as a separate DOF.  If aspects of a measurement are correlated, they should be combined
    // into a single residual
    unsigned retval(0);
    for(unsigned ires=0; ires< nResid(); ires++)
      if(activeRes(ires)) retval++;
    return retval;
  }

  template <class KTRAJ> void ResidualHit<KTRAJ>::updateWeight(MetaIterConfig const& miconfig) {
    // start with a null weight
    weight_ = Weights();
    for(unsigned ires=0; ires< nResid(); ires++) {
      if(activeRes(ires)) {
        auto const& res = residual(ires);
        // convert derivatives vector to a Nx1 matrix
        ROOT::Math::SMatrix<double,NParams(),1> dRdPM;
        dRdPM.Place_in_col(res.dRdP(),0,0);
        // convert the variance into a 1X1 matrix
        ROOT::Math::SMatrix<double, 1,1, ROOT::Math::MatRepSym<double,1>> RVarM;
        // weight by inverse variance
        double tvar = res.variance();
        RVarM(0,0) = 1.0/tvar;
        // expand these into the weight matrix
        DMAT wmat = ROOT::Math::Similarity(dRdPM,RVarM);
        // translate residual value into weight vector WRT the reference parameters
        // sign convention reflects resid = measurement - prediction
        DVEC wvec = wmat*HIT::referenceParameters().parameters() + res.dRdP()*res.value()/tvar;
        // weights are linearly additive
        weight_ += Weights(wvec,wmat);
      }
    }
    // now scale by the temp
    weight_ *= 1.0/miconfig.varianceScale();
  }
}

#endif
