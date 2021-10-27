#ifndef KinKal_MatEnv_SimpleFileFinder
#define KinKal_MatEnv_SimpleFileFinder
//
//  Simple implementation of FileFinderInterface
//
#include "KinKal/MatEnv/FileFinderInterface.hh"

namespace MatEnv {
  class SimpleFileFinder : public FileFinderInterface {
    public:
      // Specify the project and directory
      SimpleFileFinder(std::string const& project, std::string const& dir);
      SimpleFileFinder();
      virtual std::string matElmDictionaryFileName() const override { return findFile(elemlist_); }
      virtual std::string matIsoDictionaryFileName() const override { return findFile(isolist_); }
      virtual std::string matMtrDictionaryFileName() const override { return findFile(matlist_); }
      auto const& project() const { return project_; }

      // Find the specified file in the standard search path.
      virtual std::string findFile( std::string const& file ) const override { return project_ + file; }
      // default implementation
    private:
      std::string project_;
      std::string elemlist_, isolist_, matlist_;
  };
}
#endif
