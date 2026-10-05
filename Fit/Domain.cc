#include "KinKal/Fit/Domain.hh"
namespace KinKal {
  Domain::Domain(Domain const& rhs):
    range_(rhs.range()),
    bnom_(rhs.bnom()){
    }

  std::shared_ptr<Domain> Domain::clone(CloneContext& context) const{
    auto rv = std::make_shared<Domain>(*this);
    return rv;
  }
}
