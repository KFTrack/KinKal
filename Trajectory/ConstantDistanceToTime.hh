#include "KinKal/Trajectory/DistanceToTime.hh"
#pragma once

class ConstantDistanceToTime : public DistanceToTime {
    public:
        ConstantDistanceToTime(double constantSpeed, double timeOffset);
        double distance(double deltaT) override;
        double time(double distance) override;
        double speed(double distance) override;
        double inverseSpeed(double distance) override;
    private:
        double constantSpeed_;
};