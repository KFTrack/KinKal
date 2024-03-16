#ifndef KinKal_Config_hh
#define KinKal_Config_hh
//
//  This class defines the configuration of the KKTrk Fit.  It provides the external information needed by
//  the fit (BFieldMap maps, material model), the iteration schedule through which the fit should proceed,
//  and the convergence criteria for that.
//
//
#include "KinKal/General/Vectors.hh"
#include "KinKal/Fit/MetaIterConfig.hh"
#include <vector>
#include <memory>
#include <algorithm>
#include <ostream>
#include <istream>

namespace KinKal {
  struct Config {
    enum printLevel{none=0,minimal, basic, complete, detailed, extreme};
    using Schedule = std::vector<MetaIterConfig>;
    explicit Config() {}
    explicit Config(Schedule const& schedule) : Config() { schedule_ = schedule; }
    Schedule const& schedule() const { return schedule_; }
    // algebraic iteration parameters
    unsigned maxniter_ = 10; // maximum number of algebraic iterations for this config
    double dwt_ = 1.0e6; // dweighting of initial seed covariance
    double convdchisq_ = 1.0e-2; // maximum change in chisquared/dof for convergence
    double divdchisq_ = 1.0e1; // minimum change in chisquared/dof for divergence
    double pdchisq_ = 1.0e6; // maximum allowed parameter change (units of chisqred) WRT previous reference
    double divgap_ = 1.0e2; // maximum average gap of trajectory before calling it diverged (mm)
    double tol_ = 1.0e-4; // tolerance on fractional momentum accuracy due to BField domain steps
    unsigned minndof_ = 5; // minimum number of DOFs to continue fit
    bool bfcorr_ = true; // whether to make BFieldMap corrections in the fit
    bool ends_ = true; // process the passive effects at each end of the track after schedule completion
    printLevel plevel_ = none; // print level
    // schedule of meta-iterations.  These will be executed sequentially until completion or failure
    Schedule schedule_;
  };
  std::ostream& operator <<(std::ostream& os, Config const& kkconfig );
}
#endif
