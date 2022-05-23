//
//  Simple implementation of ScintHit as an array of transversely segmented constant length scintilators with the light
//  propagation axis allong the z direction.
//
#include "KinKal/Detector/ScintHit.hh"
#include <stdexcept>
namespace KinKal {

  template <class KTRAJ> class SimpleScintHit : public ScintHit<KTRAJ> {
    public:
      using HIT = Hit<KTRAJ>;
      using SCINTHIT = ScintHit<KTRAJ>;
      using PCA = PiecewiseClosestApproach<KTRAJ,Line>;
      using CA = ClosestApproach<KTRAJ,Line>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;

      SimpleScintHit(PCA const& pca,
      // ScintHit interface implementations
      virtual ~SimpleScintHit(){}
    private:
  };


}
