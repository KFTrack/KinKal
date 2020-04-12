#ifndef KinKal_KKEffBase_hh
#define KinKal_KKEffBase_hh
// base class for KKEff, for untemplated functions and common content
#include "KinKal/TDir.hh"
#include <vector>
#include <array>
#include <string>
namespace KinKal {
  class KKEffBase {
    public:
      enum Status{unprocessed=-1,processed,updated,failed};
      static std::string const& statusName(Status status);
      // common properties of all effects 
      Status status(TDir tdir) const { return status_[static_cast<std::underlying_type<TDir>::type>(tdir)]; }
      KKEffBase() : status_{{unprocessed,unprocessed}} {}
      bool wasProcessed(TDir tdir) const { return status(tdir) == processed; }
      virtual ~KKEffBase(){}
    protected:
      void setStatus(TDir tdir, Status status) { status_[static_cast<std::underlying_type<TDir>::type>(tdir)] = status; }
      void updateStatus() { status_[0] = status_[1] = updated; }
      std::array<Status,2> status_; // status of processing in each direction
      static std::vector<std::string> statusNames_;

  };
}
#endif
