#include "KinKal/TRange.hh"
namespace KinKal {
  std::ostream& operator <<(std::ostream& ost, TRange const& trange) {
    if(trange.infinite())
      ost << " Infinite Range ";
    else
      ost << " Range [" << trange.low() << "," << trange.high() << "]";
    return ost;
  }
}
