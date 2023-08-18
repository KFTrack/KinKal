#include "KinKal/Examples/CaloDistanceToTime.hh"
#include <cmath>
#include <stdexcept>
#include <cstdlib>
#include <string>

CaloDistanceToTime::CaloDistanceToTime(double asymptoticSpeed, double distanceOffset) : 
  asymptoticSpeed_(asymptoticSpeed), distanceOffset_(distanceOffset), timeOffset_(sqrt(1+pow(distanceOffset/asymptoticSpeed, 2))) {}

  double CaloDistanceToTime::distance(double deltaT) {
    // Use linear taylor approximation at t=0 if invalid
    if (deltaT > timeOffset_-1) {
        return deltaT * (pow(asymptoticSpeed_, 2) + pow(distanceOffset_, 2)) / (distanceOffset_ * timeOffset_);
    }
    return distanceOffset_ - asymptoticSpeed_ * sqrt(pow(deltaT - timeOffset_, 2) - 1);
}

double CaloDistanceToTime::time(double distance) {
    //double static const calorimeterLength = 200;
    if (distance >= distanceOffset_) {
        return timeOffset_-1;
    } else if (distance <= 0) {
        return 0;
    }
    return timeOffset_ - evaluate_root(distance);
}

double CaloDistanceToTime::speed(double distance) {
    double static const speedOfLight = 299792458.0;
    double actualSpeed = 1 / inverseSpeed(distance);
    if (abs(actualSpeed) > speedOfLight) {
        return speedOfLight;
    }
    return actualSpeed;
}

double CaloDistanceToTime::inverseSpeed(double distance) {
    if (distance >= distanceOffset_) {
        return 0;
    }
    return (distanceOffset_ - distance) / (pow(asymptoticSpeed_, 2) * evaluate_root(distance));
}

double CaloDistanceToTime::evaluate_root(double distance) {
    return sqrt(1 + pow((distance - distanceOffset_) / asymptoticSpeed_, 2));
}