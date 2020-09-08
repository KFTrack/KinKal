#ifndef KinKal_TPocaData_hh
#define KinKal_TPocaData_hh
//
// data payload for POCA calculations
//
#include "KinKal/Vectors.hh"
#include <string>
#include <vector>
namespace KinKal {

  struct TPocaData {
    enum TPStat{converged=0,unconverged,oscillating,diverged,pocafailed,invalid};
    static std::string const& statusName(TPStat status);
    //accessors
    Vec4 const& particlePoca() const { return partPoca_; }
    Vec4 const& sensorPoca() const { return sensPoca_; }
    double particleToca() const { return partPoca_.T(); }
    double sensorToca() const { return sensPoca_.T(); }
    Vec3 const& particleDirection() const { return pdir_; }
    Vec3 const& sensorDirection() const { return sdir_; }
    TPStat status() const { return status_; }
    std::string const& statusName() const { return statusName(status_); }
    double doca() const { return doca_; } // DOCA signed by angular momentum
    double docaVar() const { return docavar_; } // uncertainty on doca due to particle trajectory parameter uncertainties (NOT sensory uncertainties)
    double tocaVar() const { return tocavar_; } // uncertainty on toca due to particle trajectory parameter uncertainties (NOT sensory uncertainties)
    double dirDot() const { return pdir_.Dot(sdir_); } 
    // utility functions
    Vec4 delta() const { return sensPoca_-partPoca_; } // measurement - prediction convention
    double deltaT() const { return sensPoca_.T() - partPoca_.T(); }
    bool usable() const { return status_ < diverged; }
    TPocaData() : status_(invalid), doca_(-1.0), docavar_(-1.0), tocavar_(-1.0)  {}
    TPStat status_; // status of computation
    double doca_, docavar_, tocavar_;
    Vec3 pdir_, sdir_; // particle and sensor directions at POCA, signed by time propagation
    Vec4 partPoca_, sensPoca_; //POCA for particle and sensor
    void reset() {status_ = unconverged;}
    static std::vector<std::string> statusNames_;
  };

}
#endif
