#ifndef KinKal_TRange_hh
#define KinKal_TRange_hh
// simple struct to describe a time range
namespace KinKal {
  struct TRange {
    std::array<double,2> range_; // range of times
    TRange() : range_{1.0,-1.0} {} // initialize to have infinite range
    TRange(double low, double high) : range_{low,high} {}
    // range is define [)
    bool inRange(double t) const { return (range_[0] > range_[1]) ||
      (t >= range_[0] && t < range_[1]); }
    double low() const { return range_[0]; }
    double high() const { return range_[1]; }
    double& low() { return range_[0]; }
    double& high() { return range_[1]; }
    bool overlaps(TRange const& other ) const {
      return (high() > other.low() || low() < other.high()); }
    bool contains(TRange const& other) const {
      return (low() < other.low() && high() > other.high()); }
    static constexpr double tbuff_ = 1.0e-6; // small buffer to prevent overlaps between adjacent trajs
  };
}
#endif
