#ifndef KinKal_CylBMap_hh
#define KinKal_CylBMap_hh
// "dumb" data class to do load and store Br, Bz grid values. This object will be used by CylBFieldMap interface.
#include <string>
#include<vector>

namespace KinKal {
  class CylBMap {
    public:
      using VEC = std::vector<double>;
      using MAT = std::vector<std::vector<double>>;
      CylBMap(const std::string& file);
    private:
      int m_, n_, ntot_;
      VEC r_, z_;
      MAT Br_, Bz_;
      friend class CylBFieldMap;
  };
}
#endif