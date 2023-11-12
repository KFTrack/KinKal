
#ifndef KinKal_ParticleTrajectoryInfo_hh
#define KinKal_ParticleTrajectoryInfo_hh
#include <vector>
namespace KinKal {
  struct ParticleTrajectoryInfo {
    ParticleTrajectoryInfo(){};
    Float_t time_, dperp_, dt_;
    static std::string leafnames() { return std::string("time/f:dperp/f:dt/f"); }
  };
  typedef std::vector<ParticleTrajectoryInfo> KTIV;
}
#endif
