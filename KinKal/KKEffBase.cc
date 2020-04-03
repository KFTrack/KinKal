#include "KinKal/KKEffBase.hh"
namespace KinKal {
  std::vector<std::string> KKEffBase::statusNames_ = { "Unprocessed", "Processed", "Updated", "Failed" };
  std::string const& KKEffBase::statusName(KKEffBase::Status status) {
    switch (status) {
      case KKEffBase::unprocessed: default:
	return statusNames_[0];
      case KKEffBase::processed:
	return statusNames_[1];
      case KKEffBase::updated:
	return statusNames_[2];
      case KKEffBase::failed:
	return statusNames_[3];
    }
  }
}
