#ifndef KinKal_KKEff_hh
#define KinKal_KKEff_hh
//
// Class representing a discrete effect along the kinematic Kalman filter track fit
// This is a base class for specific subclasses representing measurements, material interactions, etc.
// Templated on the trajectory class representing the particle in this fit
//
#include "KinKal/Types.hh"
#include "KinKal/TDir.hh"
#include "KinKal/PKTraj.hh"
#include "KinKal/WData.hh"
#include <array>
#include <memory>

namespace KinKal {
// base class for KKEff, for untemplated functions and content
  class KKEffBase {
    public:
      enum Status{unprocessed=0,processed,reset,failed};
      // time of this site 
      virtual double time() const = 0;
      virtual unsigned nDOF() const = 0;
      Status status(TDir tdir) const { return status_[static_cast<std::underlying_type<TDir>::type>(tdir)]; }
      bool isActive() const { return active_; }
      KKEffBase(bool active=true) : status_{{unprocessed,unprocessed}},active_(active) {}
      virtual ~KKEffBase(){}
    protected:
      void setStatus(TDir tdir, Status status) { status_[static_cast<std::underlying_type<TDir>::type>(tdir)] = status; }
      std::array<Status,2> status_; // status of processing in each direction
      bool active_; // activity flag
  };

  template<class KTRAJ> class KKEff : public KKEffBase {
    public:
      // type of the data payload used in this trajectory
      typedef typename KTRAJ::PDATA PDATA; // forward the parameter type 
      typedef WData<KTRAJ::PDATA::PDim()> WDATA; // declare the associated weights type as well 
      typedef PKTraj<KTRAJ> PKTRAJ;
      // process this site given the adjacent site
      virtual bool process(KKEff const& other,TDir tdir) = 0;
      // update this site for a new refernce trajectory.  This must be overriden, but the base class implementation is still useful
      virtual bool update(PKTRAJ const& ref) = 0;
      // Access the processed data from this site in a particular direction; uses lazy conversion between parameter and weights space
      PDATA const& params(TDir tdir) const {
	auto idir = static_cast<std::underlying_type<TDir>::type>(tdir);
	if(!params_[idir].matrixOK() && weights_[idir].matrixOK()) params_[idir].invert(weights_[idir]);
	return params_[idir];
      }
      WDATA const& weights(TDir tdir) const {
	auto idir = static_cast<std::underlying_type<TDir>::type>(tdir);
	if(!weights_[idir].matrixOK() && params_[idir].matrixOK()) weights_[idir].invert(params_[idir]);
	return weights_[idir];
      }
      KTRAJ const& referenceTraj() const { return *reftraj_; }
      virtual ~KKEff(){} 
    protected:
      KKEff() : KKEffBase(false), reftraj_(0) {}
      KKEff(KTRAJ const& reftraj) : reftraj_(&reftraj) {}
      // allow subclasses to change the reference: this happens during updating
      void setRefTraj(KTRAJ const& newref) { reftraj_ = &newref; }
      void setRefTraj(PKTRAJ const& newref) { reftraj_ = &(newref.nearestPiece(time())); }
      // payload
      mutable std::array<PDATA,2> params_; // params after processing
      mutable std::array<WDATA,2> weights_; // weights after processing
      KTRAJ const* reftraj_; // reference trajectory for this site
  };

}

#endif
