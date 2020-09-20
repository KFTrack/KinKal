#include "General/TimeRange.hh"
namespace KinKal {
  std::ostream& operator <<(std::ostream& ost, TimeRange const& trange) {
    if(trange.infinite())
      ost << " Infinite Range ";
    else
      ost << " Range [" << trange.begin() << "," << trange.end() << "]";
    return ost;
  }
}
