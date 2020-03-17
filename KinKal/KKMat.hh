#ifndef KinKal_KKMat_hh
#define KinKal_KKMat_hh
//
// Class to describe effect of a particle passing through discrete material on the fit (ie material transport)
// This effect adds no information content, just noise, and is processed in params space 
//
#include "KinKal/KKPEff.hh"
#include "KinKal/DMat.hh"
#include "KinKal/KTMIsect.hh"
#include "KinKal/TPoca.hh"
#include "KinKal/TDir.hh"
#include <iostream>
#include <stdexcept>
namespace KinKal {
  template<class KTRAJ> class KKMat : public KKPEff<KTRAJ> {
    public:
      typedef KKPEff<KTRAJ> KKPEFF;
      typedef PKTraj<KTRAJ> PKTRAJ;
      typedef TDPoca<PKTRAJ,TLine> TDPOCA;
      typedef typename KKPEFF::PDATA PDATA; // forward the typedef
      typedef typename KTRAJ::PDER PDER; // forward the typedef
      virtual double time() const override { return time_; }
      virtual bool isActive() const override { return active_; }
      virtual bool update(PKTRAJ const& ref) override;
      // update for materials associated with a hit
      bool update(TDPOCA const& tpoca);
      virtual ~KKMat(){}
    // create from material and POCA
      KKMat(DMat const& dmat, TDPOCA const& tdpoca, bool active = true);
      // create from just the material and a trajectory 
      KKMat(DMat const& dmat, KTRAJ const& ktraj, bool active = true); 
    private:
      DMat const& dmat_; // associated detector material
      double time_;
      bool active_;
  };

   template<class KTRAJ> KKMat<KTRAJ>::KKMat(DMat const& dmat, TDPOCA const& tdpoca, bool active) : dmat_(dmat), active_(active) {
     update(tdpoca);
   }
   
   template<class KTRAJ> KKMat<KTRAJ>::KKMat(DMat const& dmat, KTRAJ const& ktraj, bool active) : dmat_(dmat), active_(active) {
     update(ktraj);
   }

   template<class KTRAJ> bool KKMat<KTRAJ>::update(PKTRAJ const& ref) {
     // not currently implemented FIXME!
     time_ = ref.range().mid();
     this->setRefTraj(ref);
     KKEffBase::updateStatus();
     KKPEFF::resetCache();
     return false;
   }

   template<class KTRAJ> bool KKMat<KTRAJ>::update(TDPOCA const& tdpoca)  {
     // update base class info
     this->setRefTraj(tdpoca.ttraj0());
     KKEffBase::updateStatus();
     KKPEFF::resetCache();
     // define the time of this effect using POCA.  Add a small offset, this should be a parameter, FIXME!
     static double epsilon_(1e-3);
     time_ = tdpoca.t0() + epsilon_;
     // find and fill the individual material intersections given this poca
     std::vector<MIsect> misects;
     dmat_.intersect(tdpoca,misects);
     // translate these to material effects
     KTMIsect<PKTRAJ> ktmisect(tdpoca.ttraj0(),tdpoca.poca0().T(),misects);
     // loop over the momentum change basis directions, adding up the effects on parameters from each
     for(int idir=0;idir<=KInter::theta2; idir++) {
       auto mdir = static_cast<KInter::MDir>(idir);
       // get the derivatives of the parameters WRT material effects
       PDER pder;
       this->refTraj().momDeriv(mdir, time(), pder);
       // convert derivative vector to a Nx1 matrix
       ROOT::Math::SMatrix<double,KTRAJ::NParams(),1> dPdm;
       dPdm.Place_in_col(pder,0,0);
       // get the mean and variance of this momentum effect from the material
       // for now, do this in the time forwards direction only, in future this may be split
       double dmom, momvar;
       ktmisect.momEffect(TDir::forwards, mdir, dmom, momvar);
       // test
       if(isnan(momvar))throw std::runtime_error("Invalid momvar");
       // update the transport for this effect; first the parameters.  Note these are for forwards time propagation (ie energy loss)
       this->pdata_.parameters() += pder*dmom;
       // now the variance: this doesn't depend on time direction
       ROOT::Math::SMatrix<double, 1,1, ROOT::Math::MatRepSym<double,1> > MVar;
       MVar(0,0) = momvar;
       this->pdata_.covariance() += ROOT::Math::Similarity(dPdm,MVar);

     }
     return true;
   }

}
#endif
