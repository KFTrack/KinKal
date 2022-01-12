
#ifndef MATISOTOPEOBJ_HH
#define MATISOTOPEOBJ_HH 1

#include <iostream>
#include <vector>
#include <cmath>

struct ModeSigma{};
struct MeanRMS{};
class MoyalDist{
  public:

    MoyalDist(ModeSigma, double mpv, double xi, int max = 20):_mode(mpv), _sigma(xi),_kmax(max) {
      _mean = _mode + _sigma *  MFACTOR;
      _rms = 0.5 * (std::pow( M_PI * _sigma , 2) );
      setCoeffs(_kmax);
    }

    MoyalDist(MeanRMS, double Mean, double RMS, int max = 20):_mean(Mean),_rms(RMS),_kmax(max) {        
       //Variance of moyal = (pi * sigma)^2 /2
      _sigma = std::sqrt(2.0) * _rms / M_PI ; 
      _mode = _mean - _sigma * MFACTOR; 
      setCoeffs(_kmax);
    }

    double sampleAR() const;
    double sampleInvCDF(double rand) const;
    double getMean() const { return _mean; }
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
    static constexpr double EG = 0.577215664901532860606;    //Euler-gamma constant; replace this with numerics when we move to c++20
    static constexpr double MFACTOR = EG + M_LN2;
};

#endif
