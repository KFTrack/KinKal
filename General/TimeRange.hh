#ifndef KinKal_TimeRange_hh
#define KinKal_TimeRange_hh
// simple struct to describe a time range, defined [ )
#include <algorithm>
#include <ostream>
#include <array>
namespace KinKal {
  class TimeRange {
    public:
      static constexpr double tbuff_ = 1.0e-10; // small buffer to prevent overlaps between adjacent trajs
      TimeRange() : range_{1.0,-1.0} {} // initialize to have infinite range
      TimeRange(double begin, double end) : range_{begin,end} {}
      bool inRange(double t) const { return infinite() || (t >= range_[0] && t < range_[1]); }
      double begin() const { return range_[0]; }
      double end() const { return range_[1]; }
      double mid() const { return 0.5*(range_[0]+range_[1]); }
      double range() const { return (range_[1]-range_[0]); }
      double& begin() { return range_[0]; }
      double& end() { return range_[1]; }
      bool infinite() const { return end() + tbuff_ < begin(); }
      bool overlaps(TimeRange const& other ) const {
	return (end() > other.begin() || begin() < other.end()); }
      bool contains(TimeRange const& other) const {
	return (begin() < other.begin() && end() > other.end()); }
      // force time to be in range
      void forceRange(double& time) const { time = std::min(std::max(time,begin()),end()); }
      // test if a time is at the limit
      bool atLimit(double time) const { return time >= end() || time <= begin(); }
    private:
      std::array<double,2> range_; // range of times
  };
  std::ostream& operator <<(std::ostream& ost, TimeRange const& trange);
}
#endif
