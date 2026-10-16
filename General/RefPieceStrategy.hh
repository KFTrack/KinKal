#ifndef KinKal_RefPieceStrategy_hh
#define KinKal_RefPieceStrategy_hh
//
//  Strategy for choosing the trajectory piece an effect references when it is updated to a new piecewise trajectory.
//  DOCA: follow the closest approach in space across pieces until the solution lies on the piece used to compute it.
//  TOCA: compute the closest approach on the piece holding the hint time, without moving to another piece.
//
#include <ostream>
#include <stdexcept>
#include <string>
namespace KinKal {
  enum struct RefPieceStrategy{DOCA=0,TOCA};
  inline std::string refPieceStrategyName(RefPieceStrategy strategy) {
    return strategy == RefPieceStrategy::TOCA ? std::string("TOCA") : std::string("DOCA");
  }
  inline RefPieceStrategy refPieceStrategy(std::string const& name) {
    if(name == "DOCA") return RefPieceStrategy::DOCA;
    if(name == "TOCA") return RefPieceStrategy::TOCA;
    throw std::invalid_argument("Unknown RefPieceStrategy " + name);
  }
  inline std::ostream& operator <<(std::ostream& ost, RefPieceStrategy strategy) {
    ost << refPieceStrategyName(strategy);
    return ost;
  }
}
#endif
