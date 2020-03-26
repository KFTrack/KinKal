//--------------------------------------------------------------------------
// File and Version Information:
//
// Description:
//      
//      
//
// Environment:
//      Software developed for the BaBar Detector at the SLAC B-Factory.
//
// Author List:
//       Matthias Steinke
//       
//
// Copyright Information:
//
//------------------------------------------------------------------------

#ifndef BBRCOLLECTIONUTILS_HH
#define BBRCOLLECTIONUTILS_HH

#include <algorithm>
namespace MatEnv {

  struct PtrLess {
    template<class PtrType>
      bool operator()(PtrType ptr1, PtrType ptr2) const {
	return (*ptr1) < (*ptr2);
      }
  };   

  struct DeleteObject{

    template< class T >
      void operator()(const T* ptr) const {
	delete ptr;
      }
  };

  struct DeleteArray{

    template< class T >
      void operator()(const T array[]) const {
	delete[] array;
      }
  };

  /**
   *  Determines the offset of the first occurrence of a specified
   *  value in a container.  This is _not_ an STL-ish way to work;
   *  the use of std::find() and iterators instead of offsets is
   *  _strongly_ preferred, even for vectors.
   *
   *  This function is supplied only as a migration aid for previous users
   *  of the Rogue Wave vector classes.
   *
   *  It is valid for any container type C for whose iterators
   *  operator- is defined.  When a restriction in the STL supplied
   *  by RW for Sun is removed in the future (see below),
   *  it could be rewritten to work for all containers.
   */
  template < class C, class T >
    typename C::difference_type
    findIndex( const C& container, const T& value ) { 
      typename C::const_iterator found 
	= std::find( container.begin(), container.end(), value );

      // Were std::distance available, one would write:               FIXME
      //// return std::distance( container.begin(), found );

      // Meanwhile, this should only compile for iterators for which it is efficient:
      return found - container.begin();
    }

}
#endif
