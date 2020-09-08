#ifndef KinKal_PTPoca_hh
#define KinKal_PTPoca_hh
//
// TPoca class dealing with piecewise trajectories
//
#include "KinKal/TPoca.hh"
#include "KinKal/PKTraj.hh"
#include <ostream>

namespace KinKal {
  template<class KTRAJ, class STRAJ> class PTPoca {
    public:
      using PKTRAJ = PKTraj<KTRAJ>;
      using KTPOCA = TPoca<KTRAJ,STRAJ>;
      // the constructor is the only non-inherited function
      PTPoca(PKTRAJ const& pktraj, STRAJ const& straj, TPocaHint const& hint, double prec=1e-6);
      // copy the TPOCA interface.  This is ugly and a maintenance burden, but avoids inheritance problems
      TPocaData::TPStat status() const { return tpdata_.status(); }
      std::string const& statusName() const { return tpdata_.statusName(); }
      double doca() const { return tpdata_.doca(); }
      double docaVar() const { return tpdata_.docaVar(); }
      double tocaVar() const { return tpdata_.tocaVar(); }
      double dirDot() const { return tpdata_.dirDot(); }
      double deltaT() const { return tpdata_.deltaT(); }
      bool usable() const { return tpdata_.usable(); }
      double particleToca() const { return tpdata_.particleToca(); }
      double sensorToca() const { return tpdata_.sensorToca(); }
      Vec4 const& particlePoca() const { return tpdata_.particlePoca(); }
      Vec4 const& sensorPoca() const { return tpdata_.sensorPoca(); }
      Vec4 delta() const { return tpdata_.delta(); }
      Vec3 const& particleDirection() const { return tpdata_.particleDirection(); }
      Vec3 const& sensorDirection() const { return tpdata_.sensorDirection(); }
      TPocaData const& tpData() const { return tpdata_; }
      PKTRAJ const& particleTraj() const { return pktraj_; }
      size_t particleTrajIndex() const { return pindex_; }
      STRAJ const& sensorTraj() const { return straj_; }
      DVEC const& dDdP() const { return dDdP_; }
      DVEC const& dTdP() const { return dTdP_; }
      bool inRange() const { return particleTraj().inRange(particleToca()) && sensorTraj().inRange(sensorToca()); }
      double precision() const { return precision_; }
      void print(std::ostream& ost=std::cout,int detail=0) const;
    private:
      double precision_; // precision used to define convergence
      TPocaData tpdata_; // data payload of POCA calculation
      PKTRAJ const& pktraj_;
      STRAJ const& straj_;
      size_t pindex_; // indices to the local traj used in TPOCA calculation
      DVEC dDdP_; 
      DVEC dTdP_;
  };

  template<class KTRAJ, class STRAJ> PTPoca<KTRAJ,STRAJ>::PTPoca(PKTRAJ const& pktraj, STRAJ const& straj, TPocaHint const& hint, double prec) : precision_(prec), pktraj_(pktraj), straj_(straj) {
    // iteratively find the nearest piece, and POCA for that piece.  Start at hints if availalble, otherwise the middle
    static const unsigned maxiter=10; // don't allow infinite iteration.  This should be a parameter FIXME!
    unsigned niter=0;
    size_t oldindex= pktraj_.pieces().size();
    pindex_ = pktraj_.nearestIndex(hint.particleToca_);
    // copy over the hint: it needs to evolve
    TPocaHint phint = hint;
    // iterate until TPOCA is on the same piece
    bool inrange(true);
    do{
      KTPOCA tpoca(pktraj_.piece(pindex_),straj,phint,prec);
      // copy the state
      tpdata_ = tpoca.tpData();
      dDdP_ = tpoca.dDdP();
      dTdP_ = tpoca.dTdP();
      inrange = tpoca.inRange();
      // update the hint
      phint.particleToca_ = tpoca.particleToca();
      phint.sensorToca_ = tpoca.sensorToca();
      // update the piece (if needed)
      oldindex = pindex_;
      pindex_ = pktraj_.nearestIndex(tpoca.particlePoca().T());
    } while( pindex_ != oldindex && usable() && niter++ < maxiter);
    // overwrite the status if we oscillated on the piece
    if(tpdata_.status() == TPocaData::converged && niter >= maxiter)
      tpdata_.status_ = TPocaData::unconverged;
    // should test explicitly for piece oscillation FIXME!
    // should test if TPOCA is in range of the piece: if not, test TPOCA against the next piece
    // and select that if it's in range, or it has with the smallest DOCA. TODO!
    std::cout << "PTPOCA niter " << niter << " status " << status() << " in range " <<  inrange << std::endl;
    if(!inrange){
      auto const& piece = pktraj_.piece(pindex_);
      std::cout << "Out of range, TOCA = " << particleToca() << " range " << piece.range() 
      << " piece " << pindex_ << " nearest " << pktraj_.nearestIndex(particleToca())
      << " npieces " << pktraj_.pieces().size() << " full range " << pktraj_.range() 
      << std::endl;

    }
  }
  

  template<class KTRAJ, class STRAJ> void PTPoca<KTRAJ,STRAJ>::print(std::ostream& ost,int detail) const {
    ost << "PTPoca status " << statusName() << " Doca " << doca() << " +- " << sqrt(docaVar())
      << " dToca " << deltaT() << " +- " << sqrt(tocaVar()) << " cos(theta) " << dirDot() << " Precision " << precision() << std::endl;
    if(detail > 0)
      ost << "Particle Poca " << particlePoca() << " Sensor Poca " << sensorPoca() << std::endl;
    if(detail > 1)
      ost << "dDdP " << dDdP() << " dTdP " << dTdP() << std::endl;
    if(detail > 2){
      ost << "Particle ";
      particleTraj().print(ost,0);
      ost << "Sensor ";
      sensorTraj().print(ost,0);
    }

  }
}
#endif

