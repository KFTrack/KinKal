#ifndef BaBaR_ErrLog_HH
#define BaBaR_ErrLog_HH
//
// Fake out the BaBaR message logger to work for Mu2e.
//
// $Id:$
// $Author:$
// $Date:$
//
// Todo:
// 1) Pass through calls to Mu2e error logger.
// 2) Call LogWarn, Log Info etc depending on severity.
// 3) For severity = fatal, throw a framework exception.
//

#include <iostream>
namespace MatEnv {

  enum Severity {debugging=-1, trace=0, routine, warning, error, fatal};

  class ErrMsg{

    public:
      explicit ErrMsg(Severity severity)
	//  : _severity(severity)
      {
      }

      ~ErrMsg()
      {
      }

      operator std::ostream&() { return std::cout; }

      template< class T >
	ErrMsg &  operator<< ( T const & t )
	{
	  std::cout << t;
	  return *this;
	}

      ErrMsg &  operator<< ( std::ostream & f(std::ostream &) )
      {
	std::cout << f;
	return *this;
      }
      ErrMsg &  operator<< ( std::ios_base & f(std::ios_base &) )
      {
	std::cout << f;
	return *this;
      }
      void  operator<< ( void f(ErrMsg &) )
      {
	f(*this);
      }

      static Severity ErrLoggingLevel;
      // private:
      //  int _severity;
  };

  // Implementation of the global operator<< that allows "endmsg" to be used
  // with a plain ostream as a manipulator.  
  inline std::ostream& operator<<( std::ostream& os, void (* fp)(ErrMsg&) )
  {
    os<<std::endl;
    return os; 
  }

  inline void endmsg( ErrMsg & err )
  {
    err << std::endl;
  }


  // Free function to inquire about the severity level.
  inline bool ErrLogging( Severity sev ){
    return (sev >= ErrMsg::ErrLoggingLevel);
  }
}
#endif
