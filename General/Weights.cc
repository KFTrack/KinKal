#include "KinKal/General/Weights.hh"
#include "KinKal/General/Parameters.hh"
namespace KinKal {
  Weights::Weights(Parameters const& pdata) : fitdata_(pdata.fitData(),true) {}
  std::ostream& operator << (std::ostream& ost, Weights const& wdata) {
    wdata.print(ost,0);
    return ost;
  }
}
