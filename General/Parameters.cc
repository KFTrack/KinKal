#include "KinKal/General/Parameters.hh"
#include "KinKal/General/Weights.hh"
#include <stdexcept>
namespace KinKal {
  Parameters::Parameters(Weights const& wdata) : fitdata_(wdata.fitData(),true) {}

  double Parameters::delta(Parameters const& other) const {
    // difference of the parameter
    DVEC pdiff = parameters()-other.parameters();
    // sum the covariances
    DMAT csum = covariance() + other.covariance();
    // invert and contract
    if(!csum.Invert())throw std::runtime_error("Inversion failure");
    double retval = ROOT::Math::Similarity(pdiff,csum);
    return retval;
  }

  std::ostream& operator << (std::ostream& ost, Parameters const& pdata) {
    pdata.print(ost,0);
    return ost;
  }
}
