#ifndef KinKal_Config_hh
#define KinKal_Config_hh
//
//  This class defines the configuration of the KKTrk Fit.  It provides the external information needed by
//  the fit (BFieldMap maps, material model), the iteration schedule through which the fit should proceed,
//  and the convergence criteria for that.
//
// struct to define a single meta-iteration of the KKTrk fit.  Each meta-iteration configuration is held
// constant until the algebraic iteration implicit in the extended Kalman fit methodology converges.
//
#include <vector>
#include <memory>
#include <algorithm>
#include <any>
#include <ostream>
#include <istream>
#include "KinKal/General/Vectors.hh"

namespace KinKal {
  struct MetaIterConfig {
    double temp_; // 'temperature' to use in the simulated annealing (dimensionless, roughly equivalent to 'sigma')
    double convdchisq_; // maximum change in chisquared/dof for convergence
    double divdchisq_; // minimum change in chisquared/dof for divergence
    int miter_; // count of meta-iteration
    // payload for effects needing special updating; specific Effect subclasses can find their particular updater inside the vector
    std::vector<std::any> updaters_;
    MetaIterConfig() : temp_(0.0), convdchisq_(0.01), divdchisq_(10.0), miter_(-1) {}
    MetaIterConfig(std::istream& is) : miter_(-1) {
      is >> temp_ >> convdchisq_ >> divdchisq_ ;
    }
    double varianceScale() const { return (1.0+temp_)*(1.0+temp_); } // variance scale so that temp=0 means no additional variance
  };

  struct Config {
    enum printLevel{none=0,minimal, basic, complete, detailed, extreme};
    enum BFCorr {nocorr=0, fixed, variable, both }; // this should be moved out of this class TODO
    using Schedule =  std::vector<MetaIterConfig>;
    explicit Config(Schedule const& schedule) : Config() { schedule_ = schedule; }
    Config() : maxniter_(10), dwt_(1.0e6),  pdchi2_(1.0e4), tbuff_(0.1), tol_(0.1), minndof_(5), bfcorr_(fixed), plevel_(none) {} 
    Schedule& schedule() { return schedule_; }
    Schedule const& schedule() const { return schedule_; }
    static bool localBFieldCorrection(BFCorr corr) { return (corr == variable || corr == both); }
    bool localBFieldCorr() const { return localBFieldCorrection(bfcorr_); }
    // algebraic iteration parameters
    int maxniter_; // maximum number of algebraic iterations for this config
    double dwt_; // dweighting of initial seed covariance
    double pdchi2_; // maximum allowed parameter change (units of chisqred) WRT previous reference
    double tbuff_; // time buffer for final fit (ns)
    double tol_; // tolerance on position change in BFieldMap integration (mm)
    unsigned minndof_; // minimum number of DOFs to continue fit
    BFCorr bfcorr_; // how to make BFieldMap corrections in the fit
    printLevel plevel_; // print level
    // schedule of meta-iterations.  These will be executed sequentially until completion or failure
    Schedule schedule_; 
  };
  std::ostream& operator <<(std::ostream& os, Config const& kkconfig );
  std::ostream& operator <<(std::ostream& os, MetaIterConfig const& miconfig );
}
#endif
