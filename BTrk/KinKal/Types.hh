#ifndef KinKal_Types_hh
#define KinKal_Types_hh
// typedefs used throughout the KinKal package
#include "Math/Cartesian3D.h"
#include "Math/LorentzVector.h"
#include "Math/PxPyPzE4D.h"
namespace KinKal {

  typedef ROOT::Math::Cartexian3D<double> ThreeV;
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D> FourV; // double precision. can used for momentum+energy and spacetime vectors

}

