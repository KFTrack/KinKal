#include "KinKal/General/Chisq.hh"
namespace KinKal {
  std::ostream& operator <<(std::ostream& ost, Chisq const& chisq ) {
    ost << "Chisq value " <<  chisq.chisq()
      << " nDOF " << chisq.nDOF()
      << " prob " << chisq.probability();
    return ost;
  }
}
