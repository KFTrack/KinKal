#ifndef KinKal_DStraw_hh
#define KinKal_DStraw_hh
//
//  description of a local segment of a straw, including a
//  mixture for the straw, the gas, and the wire, allowing for
//  offset between the wire and the straw
//
#include "KinKal/StrawMat.hh"
#include "KinKal/DPiece.hh"
#include "KinKal/TPoca.hh"
#include "KinKal/StrawMat.hh"
#include "KinKal/TLine.hh"
namespace KinKal {
  template <class TT> class DStraw : public DPiece<TT> {
    public:
    // explicit constructor from geometry and materials
      DStraw(StrawMat const& smat, TLine const& tline) : smat_(smat), tline_(tline) {}
      // DPiece interface; first, for pieces associated with a hit
      virtual void findXings(TPocaBase const& poca,std::vector<MatXing>& mxings) const override {
	return strawmat_.findXings(poca,mxings);
      }
      // also for materials without hits, just from the trajectory
      virtual void findXings(TT const& ttraj, MatXingCol& mxings) const override {
	TPoca<TT,TLine> tpoca(ttraj,axis_);
	strawmat_.findXings(tpoca,mxings);
      }
    private:
      StrawMat const& strawmat_;
      TLine axis_; // central axis
  };
}
#endif
