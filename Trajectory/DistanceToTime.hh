#pragma once

class DistanceToTime {
    public:
        virtual double distance(double deltaT) = 0;
        virtual double time(double distance) = 0;
        virtual double speed(double distance) = 0;
        virtual double inverseSpeed(double distance) = 0;
};