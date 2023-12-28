#ifndef KinKal_ExtraConfig_hh
#define KinKal_ExtraConfig_hh
//
//  This class defines the configuration for extrapolating a fit beyond the measurement region
//
#include "KinKal/General/TimeDir.hh"

namespace KinKal {
  struct ExtraConfig {
    TimeDir xdir_ = TimeDir::forwards; // time direction to extend
    double tol_ = 1.0e-3; // tolerance on fractional momentum accuracy due to BField domain steps
    double maxdt_ = 1.0e2; // maximum time to extend the fit
  };
}
#endif
