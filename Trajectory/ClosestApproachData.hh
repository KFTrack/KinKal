#ifndef KinKal_ClosestApproachData_hh
#define KinKal_ClosestApproachData_hh
//
// data payload for CA calculations
//
#include "KinKal/General/Vectors.hh"
#include <string>
#include <vector>
#include <ostream>
namespace KinKal {

  struct ClosestApproachData {
    enum TPStat{converged=0,unconverged,oscillating,diverged,parallel,failed,invalid};
    static std::string const& statusName(TPStat status);
    //accessors
    VEC4 const& particlePoca() const { return partCA_; }
    VEC4 const& sensorPoca() const { return sensCA_; }
    double particleToca() const { return partCA_.T(); }
    double sensorToca() const { return sensCA_.T(); }
    VEC3 const& particleDirection() const { return pdir_; }
    VEC3 const& sensorDirection() const { return sdir_; }
    TPStat status() const { return status_; }
    std::string const& statusName() const { return statusName(status_); }
    double doca() const { return doca_; } // DOCA signed by angular momentum
    double docaVar() const { return docavar_; } // uncertainty on doca due to particle trajectory parameter uncertainties (NOT sensory uncertainties)
    double tocaVar() const { return tocavar_; } // uncertainty on toca due to particle trajectory parameter uncertainties (NOT sensory uncertainties)
    double dirDot() const { return pdir_.Dot(sdir_); }
    double lSign() const { return lsign_; } // sign of angular momentum
    // utility functions
    VEC4 delta() const { return sensCA_-partCA_; } // measurement - prediction convention
    double deltaT() const { return sensCA_.T() - partCA_.T(); }
    bool usable() const { return status_ < diverged; }
    ClosestApproachData() : status_(invalid), doca_(-1.0), docavar_(-1.0), tocavar_(-1.0), lsign_(0.0)  {}
    TPStat status_; // status of computation
    double doca_, docavar_, tocavar_, lsign_;
    VEC3 pdir_, sdir_; // particle and sensor directions at CA, signed by time propagation
    VEC4 partCA_, sensCA_; //CA for particle and sensor
    void reset() {status_ = unconverged;}
    const static std::vector<std::string> statusNames_;
  };
  std::ostream& operator << (std::ostream& ost, ClosestApproachData const& cadata);
}
#endif
