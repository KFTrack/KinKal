#include "KinKal/Detector/Residual.hh"
namespace KinKal {

  Weights Residual::weight(DVEC const& params, double varscale) const {
    if(active()){
      // convert derivatives vector to a Nx1 matrix
      ROOT::Math::SMatrix<double,NParams(),1> dRdPM;
      dRdPM.Place_in_col(dRdP(),0,0);
      // convert the variance into a 1X1 matrix
      ROOT::Math::SMatrix<double, 1,1, ROOT::Math::MatRepSym<double,1>> RVarM;
      // weight by inverse variance
      double mvar = measurementVariance()*varscale;
      RVarM(0,0) = 1.0/mvar;
      // expand these into the weight matrix
      DMAT wmat = ROOT::Math::Similarity(dRdPM,RVarM);
      // translate residual value into weight vector WRT the reference parameters
      // sign convention reflects resid = measurement - prediction
      DVEC wvec = wmat*params + dRdP()*value()/mvar;
      // weights are linearly additive
      return Weights(wvec,wmat);
    } else
      return Weights();
  }

  std::ostream& operator <<(std::ostream& ost, Residual const& res) {
    ost << " residual value " << res.value() << " variance " << res.variance() << " dRdP " << res.dRdP();
    return ost;
  }
}
