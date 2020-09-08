#include "KinKal/PData.hh"
#include "KinKal/WData.hh"
namespace KinKal {
  PData::PData(WData const& wdata) : tdata_(wdata.tData(),true) {}

  std::ostream& operator << (std::ostream& ost, PData const& pdata) {
    pdata.print(ost,0);
    return ost;
  }
}
