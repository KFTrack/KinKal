#ifndef KinKal_TimeRange_hh
#define KinKal_TimeRange_hh
// simple struct to describe a monotonic time range, defined [).
#include <algorithm>
#include <ostream>
#include <array>
#include <stdexcept>
namespace KinKal {
  class TimeRange {
    public:
      TimeRange() : range_{0.0,0.0} {} // default to a null range (matches no times)
      TimeRange(double begin, double end) : range_{begin,end} {
        if(begin > end)throw std::invalid_argument("Invalid Time Range"); }
      double begin() const { return range_[0]; }
      double end() const { return range_[1]; }
      double rbegin() const { return range_[1]; }
      double rend() const { return range_[0]; }
      double mid() const { return 0.5*(begin()+end()); }
      double range() const { return (end()-begin()); }
      bool inRange(double t) const {return t >= begin() && t < end(); }
      bool null() const { return end() == begin(); }
      bool overlaps(TimeRange const& other ) const {
        return {inRange(other.begin()) || inRange(other.end()) || contains(other) || other.contains(*this); }
      bool contains(TimeRange const& other) const {
        return (begin() <= other.begin() && end() >= other.end()); }
      // force time to be in range
      void forceRange(double& time) const { time = std::min(std::max(time,begin()),end()); }
      // augment another range
      void combine(TimeRange const& other ) {
        range_[0] = std::min(begin(),other.begin());
        range_[1] = std::max(end(),other.end());
      }
      // restrict the range to the overlap with another range. If there is no overlap
      // leave the object unchanged and return 'false'
      bool restrict(TimeRange const& other ) {
        bool retval = overlaps(other);
        if(retval){
          range_[0] = std::max(begin(),other.begin());
          range_[1] = std::min(end(),other.end());
        }
        return retval;
      }
    private:
      std::array<double,2> range_; // range of times
  };
  std::ostream& operator <<(std::ostream& ost, TimeRange const& trange);
}
#endif
