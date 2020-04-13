#ifndef KinKal_KKConfig_hh
#define KinKal_KKConfig_hh
//
//  This class defines the configuration of the KKTrk Fit.  It provides the external information needed by
//  the fit (BField maps, material model), the iteration schedule through which the fit should proceed,
//  and the convergence criteria for that.
//
// struct to define a single meta-iteration of the KKTrk fit.  Each meta-iteration configuration is held
// constant until the algebraic iteration implicit in the extended Kalman fit methodology converges.
//
#include "KinKal/BField.hh"

#include <vector>
#include <memory>
#include <algorithm>
#include <any>
#include <ostream>

namespace KinKal {
  struct MConfig {
    float temp_; // 'temperature' to use in the simulated annealing (dimensionless, roughly equivalent to 'sigma')
    bool processmat_; // model material effect
    bool processbfield_; // model magnetic field inhomogeneity effects
    bool updatehits_; // update the internal state of the hits (activity, LR ambiguity) 
    int miter_; // count of meta-iteration
    // payload for hit updating; specific hit classes should find their particular payload inside the vector
    std::vector<std::any> hitupdateparams_;
    MConfig() : temp_(0.0), processmat_(false), processbfield_(false), updatehits_(false), miter_(-1) {}
  };

  struct KKConfig {
    typedef std::vector<MConfig> MCONFIGCOL;
    KKConfig(BField const& bfield,std::vector<MConfig>const& schedule) : bfield_(bfield),  maxniter_(10), dwt_(1.0e6), convdchisq_(0.1), divdchisq_(100.0), oscdchisq_(1.0), tbuff_(0.5), dtol_(0.1), ptol_(0.1), minndof_(5), schedule_(schedule) {} 
    KKConfig(BField const& bfield) : bfield_(bfield),  maxniter_(10), dwt_(1.0e6), convdchisq_(0.1), divdchisq_(100.0), oscdchisq_(1.0), tbuff_(0.5), dtol_(0.1), ptol_(0.1), minndof_(5){} 
    BField const& bfield() const { return bfield_; }
    MCONFIGCOL const& schedule() const { return schedule_; }
    BField const& bfield_;
    // algebraic iteration parameters
    int maxniter_; // maximum number of algebraic iterations for this config
    float dwt_; // dweighting of initial seed covariance
    float convdchisq_; // maximum change in chisquared for convergence
    float divdchisq_; // minimum change in chisquared for divergence
    float oscdchisq_; // maximum change in chisquared for oscillation
    float tbuff_; // time buffer for final fit (ns)
    float dtol_; // tolerance on direction change in BField integration (dimensionless)
    float ptol_; // tolerance on position change in BField integration (mm)
    unsigned minndof_; // minimum number of DOFs to continue fit
    // schedule of meta-iterations.  These will be executed sequentially until completion or failure
    MCONFIGCOL schedule_; 
  };
  std::ostream& operator <<(std::ostream& os, KKConfig kkconfig );
}
#endif
