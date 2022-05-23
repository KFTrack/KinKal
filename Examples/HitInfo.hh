#ifndef KinKal_HitInfo_hh
#define KinKal_HitInfo_hh
#include <vector>
#include "Math/Vector3D.h"
namespace KinKal {
  struct HitInfo {
    enum htype{straw=0, scint, parcon,unknown};
    HitInfo():
      active_(-1), type_(-1), state_(-1), ndof_(-1), id_(-1),
      time_(0.0), tresid_(0.0), tresidvar_(0.0), tresidpull_(0.0),
      dresid_(0.0), dresidvar_(0.0), dresidpull_(0.0), chisq_(0.0), prob_(0.0),
      doca_(0.0), deltat_(0.0), docavar_(0.0), tocavar_(0.0), t0_(0.0) {}

    Int_t active_, type_, state_, ndof_, id_;
    Float_t time_;
    Float_t tresid_, tresidvar_, tresidpull_, dresid_, dresidvar_, dresidpull_;
    Float_t chisq_, prob_;
    Float_t doca_, deltat_, docavar_, tocavar_, t0_;
    ROOT::Math::XYZVectorF pos_;
  };
  typedef std::vector<HitInfo> KKHIV;
}
#endif
