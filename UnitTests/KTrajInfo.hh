
#ifndef KinKal_KTrajInfo_hh
#define KinKal_KTrajInfo_hh
#include <vector>
namespace KinKal {
  struct KTrajInfo {
    KTrajInfo(){};
    ~KTrajInfo(){};
    Float_t time_, dperp_, dt_;
    static std::string leafnames() { return std::string("time/f:dperp/f:dt/f"); }
  };
  typedef std::vector<KTrajInfo> KTIV;
}
#endif
