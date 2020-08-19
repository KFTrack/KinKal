
#ifndef KinKal_KKBFieldInfo_hh
#define KinKal_KKBFieldInfo_hh
#include <vector>
namespace KinKal {
  struct KKBFieldInfo {
    KKBFieldInfo(){};
    ~KKBFieldInfo(){};
    Int_t active_;
    Float_t time_, dp_, range_;
    static std::string leafnames() { return std::string("active/i:time/f:dp/f:range/f"); }
  };
  typedef std::vector<KKBFieldInfo> KKBFIV;
}
#endif
