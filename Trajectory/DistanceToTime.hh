//
//  Define generic distance-to-time interface used to describe signal propagatation along a sensor trajectory
//  Original author: Ryan Cheng
//
#ifndef KinKal_DistanceToTime_hh
#define KinKal_DistanceToTime_hh

class DistanceToTime {
  public:
    virtual ~DistanceToTime(){}
    virtual double distance(double deltaT) = 0;
    virtual double time(double distance) = 0;
    virtual double speed(double distance) = 0;
    virtual double inverseSpeed(double distance) = 0;
};
#endif
