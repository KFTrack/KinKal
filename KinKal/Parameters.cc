#include "KinKal/Parameters.hh"
#include "KinKal/Weights.hh"
namespace KinKal {
  Parameters::Parameters(Weights const& wdata) : fitdata_(wdata.fitData(),true) {}

  std::ostream& operator << (std::ostream& ost, Parameters const& pdata) {
    pdata.print(ost,0);
    return ost;
  }
}
