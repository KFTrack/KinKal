#include "KinKal/Parameters.hh"
#include "KinKal/Weights.hh"
namespace KinKal {
  Parameters::Parameters(Weights const& wdata) : tdata_(wdata.tData(),true) {}

  std::ostream& operator << (std::ostream& ost, Parameters const& pdata) {
    pdata.print(ost,0);
    return ost;
  }
}
