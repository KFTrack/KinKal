#ifndef KinKal_MetaIterConfig_hh
#define KinKal_MetaIterConfig_hh
//
// struct to define a single meta-iteration of the KKTrk fit.  Each meta-iteration configuration is held
// constant until the algebraic iteration implicit in the extended Kalman fit methodology converges.
//
#include <vector>
#include <any>
#include <ostream>
#include <typeinfo>
#include <iomanip>
namespace KinKal {
  class MetaIterConfig {
    public:
      MetaIterConfig() : temp_(0.0) {}
      MetaIterConfig(double temp) : temp_(temp) {}
      // add updater
      void addUpdater(std::any const& updater) { updaters_.push_back(updater); }
      // accessors
      double temperature() const { return temp_; } // dimensionless parameter interpreted WRT a maximum variance
      size_t nUpdaters() const { return updaters_.size(); }
      // find a particular updater: note that at most 1 of a type is allowed for a given meta-iteration
      template<class UPDATER> const UPDATER* findUpdater() const;
    private:
      double temp_; // 'temperature' to use in the simulated annealing (dimensionless, roughly equivalent to 'sigma')
      // payload for effects needing special updating; specific Effect subclasses can find their particular updater inside the vector
      std::vector<std::any> updaters_;
  };
  template<class UPDATER> const UPDATER* MetaIterConfig::findUpdater() const {
    const UPDATER* retval(nullptr);
    for(auto const& updater : updaters_){
      auto const* myupdater = std::any_cast<UPDATER>(&updater);
      if(myupdater != 0){
        char line[100];
        if(retval != nullptr){
          snprintf(line,100,"Multiple Updaters of type %s found",typeid(UPDATER).name());
          throw std::invalid_argument(line);
        }
        retval = myupdater;
      }
    }
    return retval;
  }
  std::ostream& operator <<(std::ostream& os, MetaIterConfig const& miconfig );
}
#endif
