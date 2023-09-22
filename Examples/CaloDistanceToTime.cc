#include "KinKal/Examples/CaloDistanceToTime.hh"
#include <cmath>
#include <stdexcept>
#include <cstdlib>
#include <string>

CaloDistanceToTime::CaloDistanceToTime(double asymptoticSpeed, double distanceOffset, double flatSlope) : 
  asymptoticSpeed_(asymptoticSpeed), distanceOffset_(distanceOffset), 
  timeOffset_(sqrt(1+pow(distanceOffset/asymptoticSpeed, 2))), flatSlope_(flatSlope), 
  slopeRoot_(sqrt(1-pow(flatSlope * asymptoticSpeed, 2))) {}

  double CaloDistanceToTime::distance(double deltaT) {
    if (deltaT > timeOffset_ - 1/slopeRoot_) {
        return ((deltaT - timeOffset_ + 1/slopeRoot_) / flatSlope_) + distanceOffset_ - flatSlope_*pow(asymptoticSpeed_, 2)/slopeRoot_;
    }
    return distanceOffset_ - asymptoticSpeed_ * sqrt(pow(deltaT - timeOffset_, 2) - 1);
}

double CaloDistanceToTime::time(double distance) {
    //double static const calorimeterLength = 200;
    double cutoff = distanceOffset_ - flatSlope_*pow(asymptoticSpeed_, 2) / slopeRoot_;
    if (distance >= cutoff) {
        return (distance - cutoff) * flatSlope_ + timeOffset_-1/slopeRoot_;
    }
    return timeOffset_ - evaluate_root(distance);
}

double CaloDistanceToTime::speed(double distance) {
    // double static const speedOfLight = 299.792458; // mm/ns
    /*double actualSpeed = 1 / inverseSpeed(distance);
    if (abs(actualSpeed) > speedOfLight) {
        return speedOfLight;
    } */
    double invSpeed = inverseSpeed(distance);
    /*if (invSpeed < (1/speedOfLight)) {
        return speedOfLight;
    }*/
    return 1/invSpeed;
}

double CaloDistanceToTime::inverseSpeed(double distance) {
    if (distance >= distanceOffset_ - flatSlope_*pow(asymptoticSpeed_, 2) / slopeRoot_) {
        return flatSlope_;
    }
    return (distanceOffset_ - distance) / (pow(asymptoticSpeed_, 2) * evaluate_root(distance));
}

double CaloDistanceToTime::evaluate_root(double distance) {
    return sqrt(1 + pow((distance - distanceOffset_) / asymptoticSpeed_, 2));
}