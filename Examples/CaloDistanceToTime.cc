#include "KinKal/Examples/CaloDistanceToTime.hh"
#include <cmath>
#include <stdexcept>
#include <cstdlib>
#include <string>

CaloDistanceToTime::CaloDistanceToTime(double asymptoticSpeed, double distanceOffset) : 
  asymptoticSpeed_(asymptoticSpeed), distanceOffset_(distanceOffset), timeOffset_(sqrt(1+pow(distanceOffset/asymptoticSpeed, 2))) {}

  double CaloDistanceToTime::distance(double deltaT) {
    if (deltaT > timeOffset_-1) {
        throw std::invalid_argument("deltaT out of range with value: " + std::to_string(deltaT));
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
    double invSpeed = inverseSpeed(distance);

    if (abs(invSpeed) < 1/speedOfLight) {
        return speedOfLight;
    }
    return 1/invSpeed;
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