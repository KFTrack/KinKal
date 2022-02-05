#include "KinKal/General/TimeRange.hh"
namespace KinKal {
  std::ostream& operator <<(std::ostream& ost, TimeRange const& trange) {
    ost << " Range [" << trange.begin() << "," << trange.end() << "]";
    return ost;
  }
}
