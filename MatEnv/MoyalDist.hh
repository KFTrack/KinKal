
#ifndef MATISOTOPEOBJ_HH
#define MATISOTOPEOBJ_HH 1

#include <iostream>
#include <vector>
#include <cmath>
const double PI = 3.141592653589793238463;    //value of pi
const double EG = 0.577215664901532860606;    //Euler-gamma constant 

class MoyalDist{
    public:
        MoyalDist(double Mean, double RMS, int max = 20):_mean(Mean),_rms(RMS),_kmax(max) {
            
            //Variance of moyal = (pi * sigma)^2 /2
            _sigma = std::sqrt(2.0) * _rms / PI ; 
            //Mean = mode + sigma (( EG + std::log(2)))
            _mode = _mean - _sigma * ( EG + std::log(2)); 
            setCoeffs(_kmax);
        }
        
        double sampleAR() const;
        double sampleInvCDF(double rand) const;
        double getMean() const { return _mean; }
        double getMode() const { return _mode; }
        double getSigma()const { return _sigma; } 
        double getRMS() const { return _rms; }

    private:
        double _mode;  //For unimodal distribution most probable value and mode are the same thing 
        double _mean;  // mean of the distribution
        double _sigma; // sigma of the distribution
        double _rms;   // rms of the distribution
        int _kmax;     // Number of terms to keep in the series expansion of the inverse of cumulative distribution function 
        std::vector<double> _coeff;  // Stores the coeffs needed to calculate InvCDF for a certain _kmax
        void setCoeffs(int kmax);  // Sets the values of coeffs vector.
};

#endif