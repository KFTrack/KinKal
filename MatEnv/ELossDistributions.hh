#ifndef KinKal_MatEnv_ELossDistributions_HH
#define KinKal_MatEnv_ELossDistributions_HH

#include <iostream>
#include <vector>
#include <cmath>

namespace KinKal {
  class MoyalDist{
    public:
      struct ModeSigma{
        ModeSigma(double mpv, double xi) : mpv_(mpv), xi_(xi) {}
        double mpv_, xi_;
      };
      struct MeanRMS{
        MeanRMS(double mean, double rms) : Mean_(mean), RMS_(rms) {}
        double Mean_; double RMS_;
      };

      MoyalDist(ModeSigma const& modesigma, int max = 20):_mode(modesigma.mpv_), _sigma(modesigma.xi_),_kmax(max) {
        _mean = _mode + _sigma *  MFACTOR;
        _rms = 0.5 * (std::pow( M_PI * _sigma , 2) );
        setCoeffs(_kmax);
      }

      MoyalDist(MeanRMS const& meanrms, int max = 20):_mean(meanrms.Mean_),_rms(meanrms.RMS_),_kmax(max) {
        //Variance of moyal = (pi * sigma)^2 /2
        _sigma = std::sqrt(2.0) * _rms / M_PI ;
        _mode = _mean - _sigma * MFACTOR;
        setCoeffs(_kmax);
      }

      double sampleAR() const;
      double sample(double rand) const; // input is a random number [0,1]
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

  /// Putting in a class for Bremsstrahlung loss
  class BremssLoss
  {
  
    public:
      double sampleSTDGamma(double energy, double radthickenss) const; //standard c++ library for gamma distribution
      double sampleSSPGamma(double energy, double radthickenss) const; //implementation of gamma distribution for small shape parameter
  };
  
  
}
#endif
