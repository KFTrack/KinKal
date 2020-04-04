#ifndef KinKal_TPocaBase_hh
#define KinKal_TPocaBase_hh
//
// data payload for POCA calculations
//
#include "KinKal/Vectors.hh"
#include <string>
#include <vector>
namespace KinKal {
  class TPocaBase {
    public:
      enum TPStat{converged=0,unconverged,pocafailed,derivfailed,invalid,unknown};
      static std::string const& statusName(TPStat status);
      //accessors
      Vec4 const& particlePoca() const { return partPoca_; }
      Vec4 const& sensorPoca() const { return sensPoca_; }
      double particleToca() const { return partPoca_.T(); }
      double sensorToca() const { return sensPoca_.T(); }
      TPStat status() const { return status_; }
      std::string const& statusName() const { return statusName(status_); }
      double doca() const { return doca_; } // DOCA signed by angular momentum
      double docaVar() const { return docavar_; } // uncertainty on doca due to particle trajectory parameter uncertainties (NOT sensory uncertainties)
      double tocaVar() const { return tocavar_; } // uncertainty on toca due to particle trajectory parameter uncertainties (NOT sensory uncertainties)
      double dirDot() const { return ddot_; } // cosine of angle between traj directions at POCA
      double precision() const { return precision_; }
      // utility functions
      void delta(Vec4& ds) const { ds = sensPoca_-partPoca_; } // measurement - prediction convention
      double deltaT() const { return sensPoca_.T() - partPoca_.T(); }
      void delta(Vec3& ds) const { ds = sensPoca_.Vect()-partPoca_.Vect(); }
      bool usable() const { return status_ != pocafailed && status_ != unknown; }
    protected:
      TPStat status_; // status of computation
      double doca_, docavar_, tocavar_;
      double ddot_;
      double precision_; // precision used to define convergence
      Vec4 partPoca_, sensPoca_; //POCA for particle and sensor
      // default precision = 10 um on DOCA
      TPocaBase(double precision=0.01) : status_(invalid), doca_(-1.0), docavar_(-1.0), tocavar_(-1.0), ddot_(-1.0), precision_(precision)  {}
      void reset() {status_ = unknown;}
    private:
      static std::vector<std::string> statusNames_;
  };


}
#endif
