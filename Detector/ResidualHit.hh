#ifndef KinKal_ResidualHit_hh
#define KinKal_ResidualHit_hh
//
//  class representing a hit based on a set of uncorrelated but related (same sensor) measurements, expressed as residuals.
//  Each residual esimates the 1-dimenaional tension between the measurement and the value predicted by a reference trajectory
//
#include "KinKal/Detector/Hit.hh"
#include "KinKal/Detector/Residual.hh"
namespace KinKal {

  template <class KTRAJ> class ResidualHit : public Hit<KTRAJ> {
    public:
      // override of some Hit interface.  Subclasses must still implement update and material methods
      using HIT = Hit<KTRAJ>;
      bool active() const override { return nDOF() > 0; }
      Chisq chisq() const override;
      Chisq chisq(Parameters const& params) const override;
      // ResidualHit specific interface.
      unsigned nDOF() const;
      // describe residuals associated with this hit
      virtual unsigned nResid() const = 0;
      // individual residuals may be active or inactive
      virtual bool activeRes(unsigned ires) const = 0;
      // reference residuals for this hit.  iDOF indexs the measurement and is hit-specific, outside the range will throw
      virtual Residual const& residual(unsigned ires) const = 0;
      // residuals corrected to refer to the given set of parameters (1st-order)
      Residual residual(Parameters const& params, unsigned ires) const;
    protected:
      // allow subclasses to set the weight
      void setWeight();
 };

  template <class KTRAJ> Residual ResidualHit<KTRAJ>::residual(Parameters const& pdata,unsigned ires) const {
    auto const& resid = residual(ires);
    // compute the difference between these parameters and the reference parameters
    DVEC dpvec = pdata.parameters() - HIT::refparams_.parameters();
    // project the parameter differnce to residual space and 'correct' the reference residual to be WRT these parameters
    double uresid = resid.value() - ROOT::Math::Dot(dpvec,resid.dRdP());
    return Residual(uresid,resid.variance(),resid.dRdP());
  }

  template <class KTRAJ> Chisq ResidualHit<KTRAJ>::chisq() const {
    double chisq(0.0);
    unsigned ndof(0);
    for(unsigned ires=0; ires< nResid(); ires++) {
      if(activeRes(ires)) {
        ndof++;
        // find the reference residual
        auto const& res = residual(ires);
        // project the parameter covariance into a residual space variance
        double rvar = ROOT::Math::Similarity(res.dRdP(),HIT::refparams_.covariance());
        // add the measurement variance
        rvar +=  res.variance();
        // add chisq for this DOF
        chisq += (res.value()*res.value())/rvar;
      }
    }
    return Chisq(chisq, ndof);
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

  template <class KTRAJ> void ResidualHit<KTRAJ>::setWeight() {
    // start with a null weight
    Weights weight;
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
        DVEC wvec = wmat*HIT::refparams_.parameters() + res.dRdP()*res.value()/tvar;
        // weights are linearly additive
        weight += Weights(wvec,wmat);
      }
    }
    HIT::weight_ = weight;
  }

}

#endif
