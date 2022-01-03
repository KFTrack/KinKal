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
    enum State { inactive=-2, left=-1, null=0, right=1};  // state description
    State state_; // left-right ambiguity
    bool useDrift() const { return state_ == left || state_ == right; }
    bool active() const { return state_ != inactive; }
    double lrSign() const {
      switch (state_) {
        case left:
          return -1.0;
        case right:
          return 1.0;
        default:
          return 0.0;
      }
    }
    double maxndoca_; // effective maximum DOCA value for null ambiguity hits
    WireHitState(State state, double maxndoca) : state_(state), maxndoca_(maxndoca) {}
    WireHitState() : WireHitState(inactive,0.0) {}
  };
}
#endif
