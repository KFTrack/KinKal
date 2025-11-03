#include "KinKal/General/TimeRange.hh"
namespace KinKal {

  bool TimeRange::restrict(TimeRange const& other ) {
    bool retval = overlaps(other);
    if(retval){
      range_[0] = std::max(begin(),other.begin());
      range_[1] = std::min(end(),other.end());
    }
    return retval;
  }

  void TimeRange::extend(double time, TimeDir tdir) {
    // require the resulting range to be >0
    if(tdir == TimeDir::forwards){
      if(time > range_[0])range_[1] = time;
    } else {
      if(time < range_[1])range_[0] = time;
    }
  }

  std::ostream& operator <<(std::ostream& ost, TimeRange const& trange) {
    ost << " Range [" << trange.begin() << "," << trange.end() << "]";
    return ost;
  }
}
