#ifndef KinKal_TPOCABase_hh
#define KinKal_TPOCABase_hh
//
// Untemplated base class for POCA calculations
//
#include "BTrk/KinKal/Types.hh"
#include "BTrk/KinKal/TTraj.hh"
#include <string>
#include <vector>
namespace KinKal {
  class TPOCABase {
    public:
      enum TPStat{converged=0,unconverged,outsiderange,pocafailed,derivfailed,unknown};
      static std::string const& statusName(TPStat status);
      //accessors
      Vec4 const& poca(size_t itraj) const { return poca_[itraj]; }
      Vec4 const& poca0() const { return poca_[0]; }
      Vec4 const& poca1() const { return poca_[1]; }
      TPStat status() const { return status_; }
      std::string const& statusName() const { return statusName(status_); }
      double doca() const { return doca_; }
      double precision() const { return precision_; }
      // utility functions
      void delta(Vec4& ds) const { ds = poca_[1]-poca_[0]; }
      double dt() const { return poca_[1].T() - poca_[0].T(); }
      void doca(Vec3& ds) const { Vec3 p0(poca_[0]); Vec3 p1(poca_[1]); ds = p1-p0; }
      bool usable() const { return status_ != pocafailed && status_ != unknown; }
      virtual TTraj const& ttraj0() const { return *ttraj_[0]; }
      virtual TTraj const& ttraj1() const { return *ttraj_[1]; }
      // actions
      void testRange() {
	if(usable() && ((!ttraj_[0]->inRange(poca_[0].T())) || (!ttraj_[1]->inRange(poca_[1].T()))) )
	  status_ = outsiderange;
      }
    protected:
      TPStat status_; // status of computation
      double doca_; // geometric distance of closest approach (cached)
      double precision_; // precision used to define convergence
      std::array<Vec4,2> poca_; // spacetime points at GEOMETRIC closest approach
      std::array<const TTraj*,2> ttraj_; // base class pointers to trajs
      // override default constructor:
      // default precision = 10 um on DOCA
      TPOCABase(TTraj const& traj0,TTraj const& traj1, double precision=0.01) : status_(unknown), doca_(-1.0), precision_(precision),
      ttraj_{&traj0,&traj1} {}
      void reset() {status_ = unknown;}
    private:
      static std::vector<std::string> statusNames_;
  };
}
#endif
