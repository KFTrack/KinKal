#ifndef KinKal_SensorLine_hh
#define KinKal_SensorLine_hh
//
//  Class describing linear signal propagation in a sensor with an arbitrary distance-to-time relationship
//
#include <memory>
#include "KinKal/General/Vectors.hh"
#include "KinKal/Trajectory/DistanceToTime.hh"
#include "KinKal/Geometry/LineSegment.hh"

namespace KinKal {
  class SensorLine {
    public:
      // construct from the measurement position and time and signal propagation velocity (mm/ns).  Note the velocity points TOWARDS the sensor
      SensorLine(VEC4 const& mpos, VEC3 const& svel, double length);
      SensorLine(VEC3 const& mpos, double mtime, VEC3 const& svel, double length);
      // construct from measurement position and time, and opposite (far) end of the sensor.  Signals propagate from the far end to the measurement position with the given speed
      SensorLine(VEC3 const& mpos, VEC3 const& endpos, double mtime, double speed );
      SensorLine(VEC3 const& mpos, double mtime, VEC3 const& svel, double length, std::shared_ptr<DistanceToTime> d2t);
      // accessors
      auto const& measurementPosition() const { return lineseg_.end(); } // measurement is at the end (latest time)
      double measurementTime() const { return mtime_; }
      double& measurementTime() { return mtime_; } // Fit updates need to refine this
      // early and late positions (according to signal propagation direction)
      auto const& start() const { return lineseg_.start(); }
      auto const& end() const { return lineseg_.end(); }
      //
      auto const& line() const { return lineseg_; }
      double speed(double time) const { return d2t_->speed(d2t_->distance(mtime_ - time)); } // D2T is relative to measurement end
      double length() const { return lineseg_.length(); }
      VEC3 const& direction() const { return lineseg_.direction(); }
      // time to a point
      double timeTo(VEC3 const& point) const;
      // geometric accessors
      VEC3 position3(double time) const;
      VEC4 position4(double time) const;
      VEC3 velocity(double time) const; // signal velocity
      VEC3 const& direction(double time) const { return lineseg_.direction(); }
      void print(std::ostream& ost, int detail) const;
      double timeAtMidpoint() const { return mtime_ - d2t_->time(0.5*length()); }
    private:
      double mtime_; // measurement time, at the far end of the directed line segment
      std::shared_ptr<DistanceToTime> d2t_; // distance to time relationship along the sensor
      LineSegment lineseg_; // linear representation of the sensor
  };
  std::ostream& operator <<(std::ostream& ost, SensorLine const& sline);
}
#endif
