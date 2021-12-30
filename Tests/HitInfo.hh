#ifndef KinKal_HitInfo_hh
#define KinKal_HitInfo_hh
#include <vector>
namespace KinKal {
  struct HitInfo {
    enum htype{strawtime=0,strawdistance, scint,constraint,unknown};
    HitInfo(){};
    ~HitInfo(){};
    Int_t active_, type_, state_, ndof_;
    Float_t time_, resid_, residvar_, chisq_;
    Float_t xpos_, ypos_, zpos_;
    Float_t t0_;
    static std::string leafnames() { return std::string("active/i:type/i:ambig/i:dim/i:ndof/i:time/f:resid/f:residvar/f:chisq/f:xpos/f:ypos/f:zpos/f:t0/f"); }
  };
  typedef std::vector<HitInfo> KKHIV;
}
#endif
