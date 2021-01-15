
#ifndef KinKal_BFieldInfo_hh
#define KinKal_BFieldInfo_hh
#include <vector>
namespace KinKal {
  struct BFieldInfo {
    BFieldInfo(){};
    ~BFieldInfo(){};
    Int_t active_;
    Float_t time_, dp_, range_;
    static std::string leafnames() { return std::string("active/i:time/f:dp/f:range/f"); }
  };
  typedef std::vector<BFieldInfo> KKBFIV;
}
#endif
