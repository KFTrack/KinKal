#include "KinKal/Detector/Residual.hh"
namespace KinKal {
  std::ostream& operator <<(std::ostream& ost, Residual const& res) {
    ost << " residual value " << res.value() << " variance " << res.variance() << " dRdP " << res.dRdP();
    return ost;
  }
}
