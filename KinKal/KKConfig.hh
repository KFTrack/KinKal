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
    double temp_; // 'temperature' to use in the simulated annealing (dimensionless, roughly equivalent to 'sigma')
    bool updatemat_; // update material effects
    bool updatebfcorr_; // update magnetic field inhomogeneity effects
    bool updatehits_; // update the internal state of the hits (activity, LR ambiguity) 
    double convdchisq_; // maximum change in chisquared for convergence
    double divdchisq_; // minimum change in chisquared for divergence
    double oscdchisq_; // maximum change in chisquared for oscillation
    int miter_; // count of meta-iteration
    // payload for hit updating; specific hit classes should find their particular payload inside the vector
    std::vector<std::any> hitupdateparams_;
    MConfig() : temp_(0.0), updatemat_(false), updatebfcorr_(false), updatehits_(false), convdchisq_(0.1), divdchisq_(100.0), oscdchisq_(1.0), miter_(-1) {}
  };

  struct KKConfig {
    enum printLevel{none=-1, minimal, basic, complete, detailed, extreme};
    typedef std::vector<MConfig> MCONFIGCOL;
    KKConfig(BField const& bfield,std::vector<MConfig>const& schedule) : KKConfig(bfield) { schedule_ = schedule; }
    KKConfig(BField const& bfield) : bfield_(bfield),  maxniter_(10), dwt_(1.0e6),  tbuff_(0.5), tol_(0.1), minndof_(5), addmat_(true), addbf_(true), plevel_(none) {} 
    BField const& bfield() const { return bfield_; }
    MCONFIGCOL const& schedule() const { return schedule_; }
    BField const& bfield_;
    // algebraic iteration parameters
    int maxniter_; // maximum number of algebraic iterations for this config
    double dwt_; // dweighting of initial seed covariance
    double tbuff_; // time buffer for final fit (ns)
    double tol_; // tolerance on position change in BField integration (mm)
    unsigned minndof_; // minimum number of DOFs to continue fit
    bool addmat_; // add material effects in the fit
    bool addbf_; // add BField effects in the fit
    Vec3 origin_; // nominal origin for defining BNom
    printLevel plevel_; // print level
    // schedule of meta-iterations.  These will be executed sequentially until completion or failure
    MCONFIGCOL schedule_; 
  };
  std::ostream& operator <<(std::ostream& os, KKConfig kkconfig );
  std::ostream& operator <<(std::ostream& os, MConfig mconfig );
}
#endif
