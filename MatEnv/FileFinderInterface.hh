#ifndef KinKal_MatEnv_FileFinderInterface
#define KinKal_MatEnv_FileFinderInterface
//
// Definition of an interface for an externally supplied class that
// can resolve partial filenames into an absolute path.
//
#include <string>

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

}
#endif
