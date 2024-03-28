#ifndef KinKal_Vectors_hh
#define KinKal_Vectors_hh
// vetor typedefs used throughout the KinKal package
#include "Math/Cartesian3D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/Vector2D.h"
#include "Math/SVector.h"
#include "Math/SMatrix.h"
#include <stdexcept>
namespace KinKal {
  constexpr size_t NParams() { return 6; } // kinematic fit parameter space and phase space dimension
  constexpr size_t NDim() { return 3; } // number of spatial dimensions (in our universe)
// Physical vectors (space + spacetime) in GenVector format
  using VEC3 = ROOT::Math::XYZVector; // cartesian vector
  using CYL3 = ROOT::Math::Cylindrical3D<double>; // Cylindrical vector.  Context-dependent z axis definition
  using POL3 = ROOT::Math::Polar3D<double>; // 3D polar vector
  using VEC4 = ROOT::Math::XYZTVector; // double precision spacetime vector, 4th component = time or energy
  using MOM4 = ROOT::Math::PxPyPzMVector; // 4-momentum with 4th component = mass
  using POL2 = ROOT::Math::Polar2D<double>; // 2D polar vector.  Context-dependent z axis definition
  using VEC2 = ROOT::Math::XYVector; // 2D cartesian vector
  // Algebraic representations of spatial vectors: these are not miscible with the above
  using VMAT = ROOT::Math::SMatrix<double,NDim(),NDim(),ROOT::Math::MatRepSym<double,NDim()>>;  // matrix type for spatial vector covariance
  using DPDV = ROOT::Math::SMatrix<double,NParams(),NDim(),ROOT::Math::MatRepStd<double,NParams(),NDim()>>; // parameter derivatives WRT space dimension type
  using DVDP = ROOT::Math::SMatrix<double,NDim(),NParams(),ROOT::Math::MatRepStd<double,NDim(),NParams()>>; // space dimension derivatives WRT parameter type
  using SMAT = ROOT::Math::SMatrix<double,NDim(),NDim(),ROOT::Math::MatRepSym<double,NDim()>>;  // Spatial covariance matrix
  using SVEC3 = ROOT::Math::SVector<double,NDim()>; // spatial (algebraic) vector
  using PMAT = ROOT::Math::SMatrix<double,2,2,ROOT::Math::MatRepSym<double,2>>;  // Planar covariance matrix
  using PPMAT = ROOT::Math::SMatrix<double,2,NParams(),ROOT::Math::MatRepStd<double,2,NParams()>>; // Planar projection matrix
  using SVEC6 = ROOT::Math::SVector<double,NParams()>; // type for particle state vector payload
  using PSMAT = ROOT::Math::SMatrix<double,NParams(),NParams(),ROOT::Math::MatRepStd<double,NParams(),NParams()>>; // matrix type for translating to/from state and parameters; this is not symmetric
  using RMAT = ROOT::Math::SMatrix<double,NDim(),NDim(),ROOT::Math::MatRepStd<double,NDim(),NDim()>>; // algebraic rotation matrix
  // purely algebraic vectors
  using DVEC = ROOT::Math::SVector<double,NParams()>; // data vector for parameters and weights
  using DMAT = ROOT::Math::SMatrix<double,NParams(),NParams(),ROOT::Math::MatRepSym<double,NParams()>>;  // associated matrix

}
#endif
