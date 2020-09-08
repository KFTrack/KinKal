#include "KinKal/WData.hh"
#include "KinKal/PData.hh"
namespace KinKal {
  WData::WData(PData const& pdata) : tdata_(pdata.tData(),true) {}
  std::ostream& operator << (std::ostream& ost, WData const& wdata) {
    wdata.print(ost,0);
    return ost;
  }
}
