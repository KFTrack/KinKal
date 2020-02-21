#ifndef KinKal_KKEffect_hh
#define KinKal_KKEffect_hh
//
// Class representing a discrete effect along the kinematic Kalman filter track fit
// This is a base class for specific subclasses representing measurements, material interactions, etc.
// Templated on the trajectory class representing the particle in this fit
//
#include <array>
#include "KinKal/Types.hh"
#include "KinKal/TTraj.hh"
namespace KinKal {

  template<class PTRAJ> class KKEffect {
    public:
      enum Status{unprocessed=0,processed,reset,failed};
      // type of the data payload used in this trajectory
      typedef typename PTRAJ::PDATA PDATA; // forward the parameter type 
      typedef typename WData<PTRAJ::PDATA::PDim()> WDATA; // declare the associated weight type as well 
      // process this site given the adjacent site
      virtual bool process(KKEffect const& other,TDir tdir) = 0;
      // update this site for a new refernce trajectory.  This must be overriden, but the base class implementation is still useful
      virtual void update(PTRAJ const& ref)  = 0;
      // time of this site 
      virtual double time()  = 0;
      // Access the processed data from this site in a particular direction; uses lazy conversion between parameter and weight space
      PDATA const& params(TDir dir) const {
	if(!params_[dir].valid() && weight_[dir].valid()) params_[dir].invert(weight_[dir]);
	return params_[dir];
      }
      WDATA const& weight(TDir dir) const {
	if(!weight_[dir].valid() && params_[dir].valid()) weight_[dir].inverts(params_[dir]);
	return weight_[dir];
      }
      Status status(TDir dir) const { return status_[dir]; }
      PTRAJ const& referenceTraj() const { return reftraj_; }
      bool isActive() const { return active_; }
      virtual ~KKEffect(){} 
    protected:
      KKEffect(PTRAJ const& reftraj) : reftraj_(reftraj), status_(unprocessed), active_(true) {}
      std::array<Status,2> status_; // status of processing in each direction
      std::array<PDATA,2> params_; // params after processing
      std::array<WDATA,2> weight_; // weight after processing
      PTRAJ const& reftraj_; // reference trajectory for this site
      bool active_; // activity flag
  };

  void template<> KKEffect<PTRAJ>::Update(TTraj const& ref) { reftraj_ = ref; }

}

#endif
