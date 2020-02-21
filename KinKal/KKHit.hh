#ifndef KinKal_KKHit_hh
#define KinKal_KKHit_hh
//
//  class to use information from a hit in the Kinematic fit.
//  Used as part of the kinematic Kalman fit
//
#include "KinKal/KKWeight.hh"
#include "KinKal/TrkHit.hh"
#include "KinKal/TPOCA.hh"

namespace KinKal {
  class TTraj;
  template <class PTRAJ> class KKHit : public KKWeight<PTRAJ> {
    public:
      virtual unsigned nCons() const override { return NDIM; }
      TrkHit const& hit() const { return trkhit_; }
      virtual void update(PTRAJ const& ref)  override { KKEffect<PTRAJ>::update(ref);  trkhit_.update(ref); }
      virtual unsigned nDOF() const override { return trkhit_.nDOF(); }
      virtual ~KKHit(){} 
    private:
      TrkHit const& trkhit_; // hit used for this constraint
      TPOCA<KTRAJ,TLine> tpoca_; // POCA between the reference trajectory and the line returned by the TrkHit
 };
}
#endif
