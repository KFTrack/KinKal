#ifndef KinKal_KKMat_hh
#define KinKal_KKMat_hh
//
// Class to describe effect of passing through discrete material
// This effect provides information content and is processed in params space 
//
#include "KinKal/KKPEff.hh"
#include "KinKal/TDMInter.hh"
namespace KinKal {
  template<class KTRAJ> class KKMat : public KKPEff<KTRAJ> {
    public:
      typedef typename KKPEff::PDATA PDATA; // forward the typedef
      // construct from a params
      KKMat(PDATA const& params) : params_(params) {}
      virtual double time() const override { return time_; }
      virtual bool isActive() const override { return active_; }
      virtual bool update(PKTRAJ const& ref) override;
      virtual ~KKMat(){}
    // create from a MInter (material intersection)
      KKMat(DMat const& dmat, bool active = true) : dmat_(dmat), active_(active) {}
    private:
      bool active_;
      TDMInter tdminter_;
  };

  template<> bool KKMat<KTRAJ>::update(KPKTRAJ const& ref) {
  // check if the intersection of this material has changed FIXME!
    bool retval(true);
    update();
    return retval;
  }

}
#endif
