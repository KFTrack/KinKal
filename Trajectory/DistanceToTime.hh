#pragma once

class DistanceToTime {
    public:
        DistanceToTime(double timeOffset) : timeOffset_(timeOffset) {}
        virtual double distance(double deltaT) = 0;
        virtual double time(double distance) = 0;
        virtual double speed(double distance) = 0;
        virtual double inverseSpeed(double distance) = 0;
        double timeOffset_;
};