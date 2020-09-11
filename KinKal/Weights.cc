#include "KinKal/Weights.hh"
#include "KinKal/Parameters.hh"
namespace KinKal {
  Weights::Weights(Parameters const& pdata) : tdata_(pdata.tData(),true) {}
  std::ostream& operator << (std::ostream& ost, Weights const& wdata) {
    wdata.print(ost,0);
    return ost;
  }
}
