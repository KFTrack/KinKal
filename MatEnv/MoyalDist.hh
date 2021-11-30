
#ifndef MATISOTOPEOBJ_HH
#define MATISOTOPEOBJ_HH 1

#include <iostream>
#include <vector>
#include <cmath>
const double PI = 3.141592653589793238463;    //value of pi
const double EG = 0.577215664901532860606;    //Euler-gamma constant 

class MoyalDist{
    public:
        MoyalDist(double mpv, double xi, int max) {
            sigma = xi;
            mean = mpv + sigma * ( EG + std::log(2)); 
            rms = 0.5 * (std::pow( PI * sigma , 2) );

            kmax = max;
            setCoeffs(kmax);
        }
 
        MoyalDist(double mpv, double xi):kmax(20) { //default value of kmax at 20
            sigma = xi;
            mean = mpv + sigma * ( EG + std::log(2));
            rms = 0.5 * (std::pow( PI * sigma , 2) );
            setCoeffs(kmax);
        }

        void   setCoeffs(int kmax);
        double sampleAR() const;
        double sampleInvCDF(double rand) const;

        double getMean() const { return mean; }
        double getSigma()const { return sigma; } 
        double getRMS() const { return rms; }

    private:
        double mean;  // mean of the distribution
        double sigma; // sigma of the distribution
        double rms;   // rms of the distribution
        int kmax;     // Number of terms to keep in the series expansion of the inverse of cumulative distribution function 
        std::vector<double> coeff;  // Stores the coeffs needed to calculate InvCDF for a certain kmax
};

#endif