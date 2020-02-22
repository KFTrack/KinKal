#include "KinKal/TRange.hh"
namespace KinKal {
  std::ostream& operator <<(std::ostream& ost, TRange const& trange) {
    ost << " Range: " << trange.low() << ":" << trange.high();
    return ost;
  }
}
