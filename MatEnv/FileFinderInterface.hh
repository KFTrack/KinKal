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
      SimpleFileFinder() :  project_(getenv("PACKAGE_SOURCE")) {}

      virtual std::string matElmDictionaryFileName() const override { return findFile(std::string("/MatEnv/ElementsList.data")); }
      virtual std::string matIsoDictionaryFileName() const override { return findFile(std::string("/MatEnv/IsotopesList.data")); }
      virtual std::string matMtrDictionaryFileName() const override { return findFile(std::string("/MatEnv/MaterialsList.data")); }

      // Find the specified file in the standard search path.
      virtual std::string findFile( std::string const& file ) const override { return project_ + file; }
      // default implementation
    private:
      std::string project_;
  };
}
#endif // end BaBar_FileFinderInterface
