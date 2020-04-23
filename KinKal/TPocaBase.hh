#ifndef KinKal_TPocaBase_hh
#define KinKal_TPocaBase_hh
//
// data payload for POCA calculations
//
#include "KinKal/Vectors.hh"
#include <string>
#include <vector>
namespace KinKal {
  // Hint class for TPOCA calculation. TPOCA search will start at these TOCA values, if provided.  This allows to
  // disambiguate cases with multiple solutions (like looping trajectories), or to speed up calculations when an
  // approximate answer is already known.
  struct TPocaHint{
    bool particleHint_, sensorHint_; // could have info on one, both, or none
    float particleToca_, sensorToca_; // approximate values, used as starting points for cacluations
    TPocaHint() : particleHint_(false), sensorHint_(false), particleToca_(0.0), sensorToca_(0.0) {}
  };

  class TPocaBase {
    public:
      enum TPStat{converged=0,unconverged,pocafailed,derivfailed,invalid,unknown};
      static std::string const& statusName(TPStat status);
      //accessors
      Vec4 const& particlePoca() const { return partPoca_; }
      Vec4 const& sensorPoca() const { return sensPoca_; }
      float particleToca() const { return partPoca_.T(); }
      float sensorToca() const { return sensPoca_.T(); }
      TPStat status() const { return status_; }
      std::string const& statusName() const { return statusName(status_); }
      float doca() const { return doca_; } // DOCA signed by angular momentum
      float docaVar() const { return docavar_; } // uncertainty on doca due to particle trajectory parameter uncertainties (NOT sensory uncertainties)
      float tocaVar() const { return tocavar_; } // uncertainty on toca due to particle trajectory parameter uncertainties (NOT sensory uncertainties)
      float dirDot() const { return ddot_; } // cosine of angle between traj directions at POCA
      float precision() const { return precision_; }
      // utility functions
      void delta(Vec4& ds) const { ds = sensPoca_-partPoca_; } // measurement - prediction convention
      float deltaT() const { return sensPoca_.T() - partPoca_.T(); }
      void delta(Vec3& ds) const { ds = sensPoca_.Vect()-partPoca_.Vect(); }
      bool usable() const { return status_ != pocafailed && status_ != unknown; }
      TPocaBase(float precision=1e-2) : status_(invalid), doca_(-1.0), docavar_(-1.0), tocavar_(-1.0), ddot_(-1.0), precision_(precision)  {}
    protected:
      TPStat status_; // status of computation
      float doca_, docavar_, tocavar_;
      float ddot_;
      float precision_; // precision used to define convergence
      Vec4 partPoca_, sensPoca_; //POCA for particle and sensor
      void reset() {status_ = unknown;}
    private:
      static std::vector<std::string> statusNames_;
  };


}
#endif
