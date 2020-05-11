#ifndef KinKal_KKHitInfo_hh
#define KinKal_KKHitInfo_hh
#include <vector>
namespace KinKal {
  struct KKHitInfo {
    KKHitInfo(){};
    ~KKHitInfo(){};
    Float_t resid_, residvar_, chiref_, chifit_;
    static std::string leafnames() { return std::string("resid/f:residvar/f:chiref/f:chifit/f"); }
  };
  typedef std::vector<KKHitInfo> KKHIV;
}
#endif
