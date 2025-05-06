//
// specialization for LoopHelix to search for reflections
//
#include "KinKal/Trajectory/LoopHelix.hh"
#include "KinKal/Trajectory/ParticleTrajectory.hh"
namespace KinKal {
  template<> std::shared_ptr<KinKal::LoopHelix> ParticleTrajectory<KinKal::LoopHelix>::reflection(double tstart) const {
    std::shared_ptr<KinKal::LoopHelix> retval; // null by default
    // record the direction of the particle WRT the bfield at the start
    auto starttraj = this->nearestTraj(tstart);
    double startdir = starttraj->momentum3(tstart).Dot(starttraj->bnom());
    // go to the end to see if there's a possible reflection
    auto endtraj = this->backPtr();
    double tend = endtraj->range().mid();
    double enddir = endtraj->momentum3(tend).Dot(endtraj->bnom());
    if(startdir*enddir < 0.0){
      // reflection. Binary search to find the trajectory piece just after reflection
      double tdiff = 0.5*(endtraj->range().end() - starttraj->range().begin());
      double ttest = tend - tdiff;
      auto testtraj = this->nearestTraj(ttest);
      while(testtraj != endtraj){
        double testdir =  testtraj->momentum3(ttest).Dot(testtraj->bnom());
        tdiff = 0.5*std::max(tdiff, testtraj->range().range());
        if(testdir*enddir > 0){
          endtraj = testtraj;
          ttest -= tdiff;
        } else {
          ttest += tdiff;
        }
        testtraj = this->nearestTraj(ttest);
      }
      // final test
      double testdir =  endtraj->momentum3(endtraj->range().mid()).Dot(endtraj->bnom());
      if(testdir*startdir > 0.0) throw std::invalid_argument("Inconsistant reflection search");
      retval = endtraj;
    }
    return retval;
  }
}

