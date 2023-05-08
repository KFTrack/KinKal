#include "KinKal/Trajectory/CaloDistanceToTime.hh"
#include <cmath>
#include <stdexcept>
#include <cstdlib>

CaloDistanceToTime::CaloDistanceToTime(double asymptoticSpeed, double timeOffset, double distanceOffset) : 
  DistanceToTime(timeOffset), asymptoticSpeed_(asymptoticSpeed), distanceOffset_(distanceOffset) {}

double CaloDistanceToTime::distance(double deltaT) {
    if (deltaT > timeOffset_) {
        throw std::invalid_argument("deltaT out of range with value: " + std::to_string(deltaT));
    }
    return distanceOffset_ + asymptoticSpeed_ * sqrt(pow(deltaT - timeOffset_ - 1, 2) - 1);
}

double CaloDistanceToTime::time(double distance) {
    double static const calorimeterLength = 200;
    if (distance <= distanceOffset_) {
        return timeOffset_;
    } else if (distance >= calorimeterLength) {
        return timeOffset_+1 - evaluate_root(calorimeterLength);
    }
    return timeOffset_+1-evaluate_root(distance);
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
    if (distance < distanceOffset_) {
        return 0;
    }
    return (distanceOffset_-distance) / (pow(asymptoticSpeed_, 2) * evaluate_root(distance));
}

double CaloDistanceToTime::evaluate_root(double distance) {
    return sqrt(1 + pow((distance - distanceOffset_) / asymptoticSpeed_, 2));
}