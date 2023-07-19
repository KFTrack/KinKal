//
//  Trivial implementation of Distance-to-time using a constant speed
//
#include "KinKal/Trajectory/DistanceToTime.hh"
#ifndef KinKal_ConstantDistanceToTime_hh
#define KinKal_ConstantDistanceToTime_hh

class ConstantDistanceToTime : public DistanceToTime {
  public:
    ConstantDistanceToTime(double constantSpeed);
    virtual ~ConstantDistanceToTime(){}
    double distance(double deltaT) override;
    double time(double distance) override;
    double speed(double distance) override;
    double inverseSpeed(double distance) override;
  private:
    double constantSpeed_;
};
#endif
