#ifndef KinKal_KKSite_hh
#define KinKal_KKSite_hh
//
// Class representing a processing step of the kinematic Kalman filter
// This is a base class for specific subclasses representing measurements, material interactions, etc.
// Templated on the kinematic trajectory used in this fit
//
#include <array>
#include "BTrk/KinKal/Types.hh"
namespace KinKal {

  template<class KTRAJ> class KKSite {
    public:
      enum Status{unprocessed=0,processed,failed};
      // type of the data payload used in this trajectory
      typedef typename KTRAJ::TDATA TDATA; // forward the type definition
      // simple struct to simplify storage
      struct DataPayload {
	enum DataType{param=0,weight};
	TDATA data_;
	DataType type_;
      };
      // process this site given the adjacent site
      virtual bool process(KKSite const& other,TDir tdir) = 0;
      // update this site for a new refernce.  This must be overriden, but the base class implementation is still useful
      virtual void update(KTraj const& ref)  = 0;
      // Access the processed data from this site in a particular direction
      virtual DataPayload const& processedData(TDir dir) const =0;
      Status status(TDir dir) const { return status_[dir]; }
      KTRAJ const& referenceTraj() const { return reftraj_; }
      bool isActive() const { return active_; }
    protected:
      KKSite(double time) : status_{unprocessed}, time_(time), active_(true) {}
      std::array<Status,2> status_; // status of processing in each direction
      KTRAJ const& reftraj_; // reference trajectory for this site
      double time_; // time of this site
      bool active_; // activity flag
  };


  void template<> KKSite<KTRAJ>::Update(KTraj const& ref) { reftraj_ = ref; }

}

#endif
