#ifndef BaBar_FileFinderInterface
#define BaBar_FileFinderInterface

//
// Definition of an interface for an externally supplied class that
// can resolve partial filenames into an absolute path.
//

#include <string>
#include <cstdlib>

namespace MatEnv {
  class FileFinderInterface {

    public:

      // These must be supplied.
      virtual std::string matElmDictionaryFileName() const = 0;
      virtual std::string matIsoDictionaryFileName() const = 0;
      virtual std::string matMtrDictionaryFileName() const = 0;

      // Find the specified file in the standard search path.
      virtual std::string findFile( std::string const& ) const = 0;

    protected:
      FileFinderInterface(){}
      virtual ~FileFinderInterface(){}
  };

  class SimpleFileFinder : public FileFinderInterface {
    public:
      SimpleFileFinder() :  project_(std::string(getenv("KINKAL_SOURCE_DIR")) + std::string("/MatEnv/")),
      elemlist_("ElementsList.data"),
      isolist_("IsotopesList.data"),
      matlist_("MaterialsList.data")
      {}
      // allow overriding the materials
      SimpleFileFinder(std::string const& project, std::string const& dir) :  project_(std::string(getenv(project.c_str())) + dir),
      elemlist_("ElementsList.data"),
      isolist_("IsotopesList.data"),
      matlist_("MaterialsList.data")
      {}

      virtual std::string matElmDictionaryFileName() const override { return findFile(elemlist_); }
      virtual std::string matIsoDictionaryFileName() const override { return findFile(isolist_); }
      virtual std::string matMtrDictionaryFileName() const override { return findFile(matlist_); }

      // Find the specified file in the standard search path.
      virtual std::string findFile( std::string const& file ) const override { return project_ + file; }
      // default implementation
    private:
      std::string project_;
      std::string elemlist_, isolist_, matlist_;
  };
}
#endif // end BaBar_FileFinderInterface
