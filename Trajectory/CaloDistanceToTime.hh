#include "KinKal/Trajectory/DistanceToTime.hh"
#pragma once

class CaloDistanceToTime : public DistanceToTime {
    public:
        CaloDistanceToTime(double asymptoticSpeed, double distanceOffset);
        double distance(double deltaT) override;
        double time(double distance) override;
        double speed(double distance) override;
        double inverseSpeed(double distance) override;
    private:
        double evaluate_root(double distance);
        double asymptoticSpeed_;
        double distanceOffset_;
};