
#ifndef KinKal_DomainWallInfo_hh
#define KinKal_DomainWallInfo_hh
#include <vector>
namespace KinKal {
  struct DomainWallInfo {
    DomainWallInfo(){};
    Int_t active_;
    Float_t time_, range_;
    static std::string leafnames() { return std::string("active/i:time/f:range/f"); }
  };
  typedef std::vector<DomainWallInfo> KKBFIV;
}
#endif
