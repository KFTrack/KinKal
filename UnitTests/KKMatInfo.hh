
#ifndef KinKal_KKMatInfo_hh
#define KinKal_KKMatInfo_hh
#include <vector>
namespace KinKal {
  struct KKMatInfo {
    KKMatInfo(){};
    ~KKMatInfo(){};
    Int_t active_, nxing_;
    Float_t time_, dmomf_, momvar_, perpvar_;
    static std::string leafnames() { return std::string("active/i:nxing/i:time/f:dmomf/f:momvar/f:perpvar/f"); }
  };
  typedef std::vector<KKMatInfo> KKMIV;
}
#endif
