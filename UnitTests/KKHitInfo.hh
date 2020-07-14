#ifndef KinKal_KKHitInfo_hh
#define KinKal_KKHitInfo_hh
#include <vector>
namespace KinKal {
  struct KKHitInfo {
    KKHitInfo(){};
    ~KKHitInfo(){};
    Float_t resid_, residvar_, fitchi_;
    Int_t active_;
    static std::string leafnames() { return std::string("active/i:resid/f:residvar/f:fitchi/f"); }
  };
  typedef std::vector<KKHitInfo> KKHIV;
}
#endif
