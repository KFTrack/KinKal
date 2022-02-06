#ifndef KinKal_MetaIterConfig_hh
#define KinKal_MetaIterConfig_hh
//
// struct to define a single meta-iteration of the KKTrk fit.  Each meta-iteration configuration is held
// constant until the algebraic iteration implicit in the extended Kalman fit methodology converges.
//
#include <vector>
#include <any>
#include <ostream>
namespace KinKal {
  struct MetaIterConfig {
    double temp_; // 'temperature' to use in the simulated annealing (dimensionless, roughly equivalent to 'sigma')
    // payload for effects needing special updating; specific Effect subclasses can find their particular updater inside the vector
    std::vector<std::any> updaters_;
    MetaIterConfig() : temp_(0.0) {}
    MetaIterConfig(double temp) : temp_(temp) {}
    double varianceScale() const { return (1.0+temp_)*(1.0+temp_); } // variance scale so that temp=0 means no additional variance
  };
  std::ostream& operator <<(std::ostream& os, MetaIterConfig const& miconfig );
}
#endif
