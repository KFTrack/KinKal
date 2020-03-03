#ifndef KinKal_KKMat_hh
#define KinKal_KKMat_hh
//
// Class to describe effect of a particle passing through discrete material on the fit (ie material transport)
// This effect adds no information content, just noise, and is processed in params space 
//
#include "KinKal/KKPEff.hh"
#include "KinKal/KTMIsect.hh"
#include "KinKal/TDir.hh"
namespace KinKal {
  template<class KTRAJ> class KKMat : public KKPEff<KTRAJ> {
    public:
      typedef KKPEff<KTRAJ> KKPEFF;
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef typename KKPEff::PDATA PDATA; // forward the typedef
      typedef typename KTRAJ::PDER PDER; // forward the typedef
      virtual double time() const override { return ktmisect_.tinter_;}
      virtual bool isActive() const override { return active_; }
      virtual bool update(PKTRAJ const& ref) override;
      virtual ~KKMat(){}
    // create from a TDMInter (material intersection)
      KKMat(TDMInter const& ktmisect, bool active = true);
    private:
      bool active_;
      KTMIsect<PKTRAJ> ktmisect_;
  };

   template<class KTRAJ> KKMat<KTRAJ>::KKMat(KTMIsect const& ktmisect, bool active = true) :
   ktmisect_(ktmisect), active_(active) {
     update(ktimsect.kTraj());
   }

  template<class KTRAJ> bool KKMat<KTRAJ>::update(KPKTRAJ const& ref) {
  // check if the intersection of this material has changed FIXME!

  // loop over the momentum change basis directions, adding up the effects on parameters from each
    for(int idir=0;idir<=KInter::theta2) {
      auto mdir = static_cast<KInter::MDir>(mdir);
    // get the derivatives of the parameters WRT material effects
      PDER pder;
      this->refTraj().momDeriv(tdir, time(), pder);
      // convert derivative vector to a Nx1 matrix
      ROOT::Math::SMatrix<double,KTRAJ::NParams(),1> dPdm;
      dPdm.Place_in_col(pder,0,0);
      // get the mean and variance of this momentum effect from the material
      // for now, do this in the time forwards direction only, in future this may be split
      double dmom, momvar;
      ktmisect_.momEffect(TDir::forwards, mdir, dmom, momvar);
      // update the transport for this effect; first the parameters.  Note these are for forwards time propagation (ie energy loss)
      this->pdata_.parameters() += pder*dmom;
      // now the variance: this doesn't depend on time direction
      ROOT::Math::SMatrix<double, 1,1, ROOT::Math::MatRepSym<double,1> > MVar;
      MVar(0,0) = momvar;
      this->pdata_.covariance() += ROOT::Math::Similarity(dPdm,MVar);
      // ready for processing!
    }
    return true;
  }

}
#endif
