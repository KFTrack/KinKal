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
      unsigned nDOF() const override;
      Chisq chisq(Parameters const& params) const override;
      Weights const& weight() const override { return weight_; }
      // describe residuals associated with this hit
      virtual unsigned nResid() const = 0;
      // reference residuals for this hit.  ires indexs the measurement and is hit-specific, outside the range will throw
      // this is generally biased as the refefence includes the effect of this hit
      virtual Residual const& refResidual(unsigned ires) const = 0;
     // residuals corrected to refer to the given set of parameters (1st-order)
      Residual residual(Parameters const& params, unsigned ires) const;
      // unbiased residuals WRT the reference parameters; computed from the reference
      Residual residual(unsigned ires) const;
      // unbiased pull of this residual (including the uncertainty on the reference parameters)
      double pull(unsigned ires) const;
    protected:
      ResidualHit() {}
      // ResidualHit specific interface
      void updateWeight(MetaIterConfig const& config);
    private:
      Weights weight_; // weight of this hit computed from the residuals
  };

  template <class KTRAJ> Residual ResidualHit<KTRAJ>::residual(Parameters const& params,unsigned ires) const {
    auto const& resid = refResidual(ires);
    // compute the difference between these parameters and the reference parameters
    DVEC dpvec = params.parameters() - HIT::referenceParameters().parameters();
    // project the parameter differnce to residual space and 'correct' the reference residual to be WRT these parameters
    double uresid = resid.value() - ROOT::Math::Dot(dpvec,resid.dRdP());
    double pvar = ROOT::Math::Similarity(resid.dRdP(),params.covariance());
    if(pvar<0) throw std::runtime_error("Covariance projection inconsistency");
    return Residual(uresid,resid.variance(),pvar,resid.active(),resid.dRdP());

  }

  template <class KTRAJ> Residual ResidualHit<KTRAJ>::residual(unsigned ires) const {
    return residual(HIT::unbiasedParameters(),ires);
  }

  template <class KTRAJ> double ResidualHit<KTRAJ>::pull(unsigned ires) const {
    auto ures = residual(ires);
    return ures.pull();
  }

  template <class KTRAJ> Chisq ResidualHit<KTRAJ>::chisq(Parameters const& params) const {
    double chisq(0.0);
    unsigned ndof(0);
    for(unsigned ires=0; ires< nResid(); ires++) {
      auto resid = residual(params,ires);
      chisq += resid.chisq();
      ndof += resid.nDOF();
    }
    return Chisq(chisq,ndof);
  }

  template <class KTRAJ> unsigned ResidualHit<KTRAJ>::nDOF() const {
    // each residual is counted as a separate DOF.  If aspects of a measurement are correlated, they should be combined
    // into a single residual
    unsigned retval(0);
    for(unsigned ires=0; ires< nResid(); ires++)
      retval += refResidual(ires).nDOF();
    return retval;
  }

  template <class KTRAJ> void ResidualHit<KTRAJ>::updateWeight(MetaIterConfig const& miconfig) {
    // start by zeroing the weight, then augment with each residual's weight
    weight_ = Weights();
    for(unsigned ires=0; ires< nResid(); ires++) {
      auto const& resid = refResidual(ires);
      if(resid.active())weight_ += resid.weight(HIT::referenceParameters().parameters(),miconfig.varianceScale());
    }
  }
}

#endif
