#ifndef KinKal_WireHitStructs_hh
#define KinKal_WireHitStructs_hh
#include <stdexcept>
namespace KinKal {
// struct describing local drift info
  struct DriftInfo {
    DriftInfo() : tdrift_(0.0), tdriftvar_(0.0), vdrift_(0.0) {}
    double tdrift_; // drift time
    double tdriftvar_; // variance on drift time
    double vdrift_; // instantanious drift speed
  };

  // struct describing wire hit internal state
  struct WireHitState {
    enum LRAmbig { left=-1, null=0, right=1};  // ambiguity
    enum Dimension {none=-1, time=0, distance=1, both=2}; // what gets constrained
    LRAmbig lrambig_; // left-right ambiguity
    Dimension dimension_; // physical dimensions being constrained
    double nullvar_, nulldt_ ; // spatial variance and offset for null ambiguity hits
    WireHitState(LRAmbig lrambig, Dimension dim,double nvar, double ndt) : lrambig_(lrambig), dimension_(dim), nullvar_(nvar), nulldt_(ndt) {
      if(dimension_ > time && (lrambig_ != null)) throw std::invalid_argument("Inconsistant wire hit state");
    }
    WireHitState() : WireHitState(null,none,1.0,0.0) {}
  };
}
#endif
