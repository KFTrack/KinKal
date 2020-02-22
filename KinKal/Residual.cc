#include "KinKal/Residual.hh"
namespace KinKal  {
  std::ostream& operator <<(std::ostream& ost, Residual const& res) {
    ost << "Residual " << res.resid() << " Variance " << res.residVar() << " dRdD " << res.dRdD() << " dRdT " << res.dRdT();
    return ost;
  }
}
