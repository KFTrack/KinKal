#ifndef KinKal_DPiece_hh
#define KinKal_DPiece_hh
//
//  Interface for describing a physical piece of a detector.  This implements the methods
//  for finding particle trajectory crossings of the piece.  Templated
//  on the trajectory type of the fit.  Crossing methods for specific
//  detector pieces must be specialized on the trajectory type
//
#include "KinKal/MatXing.hh"
#include "KinKal/TPocaBase.hh"
#include <memory>
namespace KinKal {
  template <class TT> class DPiece {
    public:
      virtual ~DPieceInter() {} 
      virtual void findXings(TT const& ttraj, MatXingCol& mxings) const =0;
      virtual void findXings(TPocaBase const&, MatXingCol& mxings) const =0;
  };
}
#endif
