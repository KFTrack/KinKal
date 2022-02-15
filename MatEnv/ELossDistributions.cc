#include <random>
#include <cmath>
#include <stdexcept>
#include <cmath>
#include "KinKal/MatEnv/DetMaterial.hh"
#include "KinKal/MatEnv/ELossDistributions.hh"

namespace KinKal {
  void MoyalDist::setCoeffs(int kmax){
    if(_coeff.size() == (long unsigned int) kmax){
      return; //if the coefficient vector is already filled then don't do anything
    }
    if(_kmax < 1){
      throw std::invalid_argument("kmax (number of terms in the Moyal series expansion) cannot be less than one.");
    }
    _coeff.clear();
    _coeff.push_back(1); // coeff[0] = 1;
    for(int k=0; k < kmax; k++ ){
      double sumCoeff = 0;
      if( k > 0 ){
        for(int m = 0; m < k; m++){
          sumCoeff += ( _coeff.at(m) * _coeff.at( k - 1 - m) ) / ( (m + 1) * (2.0 * m + 1) );
        }
        _coeff.push_back(sumCoeff);

      }
    }
  }

  double MoyalDist::sample(double rand) const{
    if(rand < 0 || rand >1){
      throw std::invalid_argument("random number should be between 0-1");
    }

    // The Moyal distribution cdf is cdf(x; mu, sigma) = erfc(exp(- 1/2 * (x - mu) / sigma) / sqrt(2))
    // This is an implementation of erfc^-1 (z) from https://functions.wolfram.com/06.31.06.0002.01
    // kmax is the number of terms in the series to keep. Higher the kmax better is the agreement with accept/reject m
    // but high kmax also incurs high computational cost

    double t = 0.5 * std::sqrt(M_PI) * (1 - rand);
    double sum = 0;

    for(int k=0; k < _kmax; k++ ){
      sum +=  ( _coeff.at(k) / (2.* k + 1.0) ) * std::pow(t, (2.* k + 1.0));
    }

    double y = sum;
    double x = _mode - 2 * _sigma * std::log ((std::sqrt(2.) * y));
    return x;
  }


  double MoyalDist::sampleAR() const{

    // This is an implementation of accep-reject method for sampling Moyal distribution
    // recipe at http://www.stat.rice.edu/~dobelman/textfiles/DistributionsHandbook.pdf
    // Needs two random numbers.

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);

    const double hmax = 0.912;
    double h = 0;
    double Hy = 0;
    double z = 0;

    do
    {
      double rand1 = dis(gen);
      double rand2 = dis(gen);
      //std::cout << "rand1 == " << rand1 << " && " << "rand2 == " << rand2 << "\n";

      double y = M_PI * (rand1 - 0.5);
      h = rand2 * hmax;
      z = std::tan(y);
      //std::cout << "z == " << z << "\n" ;
      Hy = std::sqrt(1.0/(2.0 * M_PI)) * (1.0 + std::pow(z, 2) ) \
           * std::exp( -0.5 * ( z + std::exp(-z) ) );

    } while ( h > Hy );
    double x = _mode + _sigma * z;
    return (x) ;
  }

  double BremssLoss::sampleSTDGamma(double energy, double radthickness) const{
    // We can use the std::gamma_distribution; However it is computaionally inefficient 
    // for very small shape parameter and generates a lot of zeros. 
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::gamma_distribution<> dis(radthickness/M_LN2, 1.0);
    double Y = 0; 
    
    Y =  dis(gen);
    
    double bremssFraction = std::exp(-Y); // [(E0 - hn)/E0]
    double bremssLoss = energy * (1. - bremssFraction);

    return bremssLoss;

  }

  double BremssLoss::sampleSSPGamma(double energy, double radthickness) const{
    // This uses a method from https://arxiv.org/pdf/1302.1884.pdf
    // which is supposed to work better for small shape parameters 
    // Performs only slightly better than std::gamma_distribution 
    double alpha = radthickness/M_LN2;
    double lambda = -1. + 1./alpha;
    double omega = alpha/(exp(1.)*(1.-alpha));
    double rate = 1./(1. + omega);
    double z =0;
    double Y =0;
    double eta_alpha = 0;
    double h_alpha = 0;
    
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);

    do
        {
           double U = dis(gen);
            if (U <= rate)
            {
                z = - log(U/rate);
            } else {
                z = log(dis(gen))/lambda;
            }

            h_alpha = exp( -z - exp(-z/alpha));
            if(z>=0){
                eta_alpha = exp(-z);
            } else{
                eta_alpha = omega * lambda * exp(lambda *z);
            }

            Y = exp(-z/alpha);
        } while (h_alpha/eta_alpha < dis(gen) );
        
        
        double bremssFraction = std::exp(-Y); // [(E0 - hn)/E0]
        double bremssLoss = energy * (1. - bremssFraction);

        return bremssLoss;

  }

  double DeltaRayLoss::sampleDRL(){

    // Since the production of a high energy DR is rare and independent of each other,
    // we assume the distrinbution of number of delta rays produced above cut-off energy
    // is Poisson with a mean of _avgNumber
    std::random_device rd;
    std::mt19937 gen(rd());
    std::poisson_distribution<> pois(_avgNumber);

    int producedDR = pois(gen);

    // Sample producedDR energies from the CDF and sum them to obtain the 
    // energy loss of each primary particle
    double elossSum = 0;
    double elossRand = 0;
    double elossRandPrevious = 0;
    int counter = 0;
    std::uniform_real_distribution<double> unif(0,1);
    for(int i=0; i < producedDR; i++){

      // The following loop is for solving the equation for inverse CDF as suggested in
      // Watts Jr, J. W. Calculation of energy deposition distributions for simple geometries. No. M452. 1973.
      do 
      { 
        if(counter == 0){
          elossRandPrevious = _cutoffEnergy;

        } else {
          elossRandPrevious = elossRand;
        }
        double rand = unif(gen);
        double elossRandInv =  ( 1. /_cutoffEnergy) + (_beta*_beta/_elossMax)*std::log(_cutoffEnergy/elossRandPrevious) \
                                - rand * ( ( 1. /_cutoffEnergy) -  ( 1. /_elossMax) + (_beta*_beta/_elossMax) * std::log(_cutoffEnergy/_elossMax) );
                             
        elossRand =  1./elossRandInv;
        counter++;

      } while (abs(elossRand - elossRandPrevious)/abs(elossRand) > 0.001);

      elossSum += elossRand;
      
    }
    return elossSum;
    
  }

}
