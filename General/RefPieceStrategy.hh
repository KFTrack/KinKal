#ifndef KinKal_RefPieceStrategy_hh
#define KinKal_RefPieceStrategy_hh
//
//  Strategy for choosing the trajectory piece an effect references when it is updated to a new piecewise trajectory.
//  DOCA: follow the closest approach in space across pieces until the solution lies on the piece used to compute it,
//        moving at most MetaIterConfig::maxPieceIter() times.
//  TOCA: compute the closest approach on the piece holding the hint time, without moving to another piece
//        (PiecewiseClosestApproach with maxiter = 0). The hint then defines the reference piece, not just a starting
//        value. Effects that take the hint from their previous solution can still move one piece per update, since
//        that solution's TOCA may lie outside its piece.
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
