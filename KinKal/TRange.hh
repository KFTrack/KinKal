#ifndef KinKal_TRange_hh
#define KinKal_TRange_hh
// simple struct to describe a time range, defined [ )
#include <algorithm>
#include <ostream>
#include <array>
namespace KinKal {
  class TRange {
    public:
      static constexpr double tbuff_ = 1.0e-10; // small buffer to prevent overlaps between adjacent trajs
      TRange() : range_{1.0,-1.0} {} // initialize to have infinite range
      TRange(double low, double high) : range_{low,high} {}
      bool inRange(double t) const { return (range_[0] > range_[1]) ||
	(t >= range_[0] && t < range_[1]); }
      double low() const { return range_[0]; }
      double high() const { return range_[1]; }
      double mid() const { return 0.5*(range_[0]+range_[1]); }
      double range() const { return (range_[1]-range_[0]); }
      double& low() { return range_[0]; }
      double& high() { return range_[1]; }
      bool infinite() const { return high() < low(); }
      bool overlaps(TRange const& other ) const {
	return (high() > other.low() || low() < other.high()); }
      bool contains(TRange const& other) const {
	return (low() < other.low() && high() > other.high()); }
      // force time to be in range
      void forceRange(double& time) const { time = std::min(std::max(time,low()),high()); }
      // test if a time is at the limit
      bool atLimit(double time) const { return time >= high() || time <= low(); }
    private:
      std::array<double,2> range_; // range of times
  };
  std::ostream& operator <<(std::ostream& ost, TRange const& trange);
}
#endif
