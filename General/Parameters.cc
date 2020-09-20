#include "General/Parameters.hh"
#include "General/Weights.hh"
namespace KinKal {
  Parameters::Parameters(Weights const& wdata) : fitdata_(wdata.fitData(),true) {}

  std::ostream& operator << (std::ostream& ost, Parameters const& pdata) {
    pdata.print(ost,0);
    return ost;
  }
}
