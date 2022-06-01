
#ifndef KinKal_MaterialInfo_hh
#define KinKal_MaterialInfo_hh
#include <vector>
namespace KinKal {
  struct MaterialInfo {
    MaterialInfo(): active_(false), nxing_(-1), time_(0.0), dmomf_(0.0),
    momvar_(0.0), perpvar_(0.0), doca_(0.0), docavar_(0.0), dirdot_(0.0){};
    ~MaterialInfo(){};
    Int_t active_, nxing_;
    Float_t time_, dmomf_, momvar_, perpvar_;
    Float_t doca_, docavar_, dirdot_;
  };
  typedef std::vector<MaterialInfo> KKMIV;
}
#endif
