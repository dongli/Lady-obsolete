#ifndef __lady_commons__
#define __lady_commons__

#include "geomtk.h"

#include <list>

using geomtk::SphereDomain;
using geomtk::SphereCoord;
using geomtk::TimeLevels;
using geomtk::RLLMesh;
using geomtk::RLLMeshIndex;

namespace lady {

#define LADY_DOMAIN SphereDomain
#define LADY_SPACE_COORD SphereCoord
#define LADY_MATRIX arma::mat
#define LADY_LIST std::list

// -----------------------------------------------------------------------------
// advection only part
#ifdef ADVECTION_ONLY
#define LADY_MESH RLLMesh
#define LADY_MESH_INDEX RLLMeshIndex
#endif

}

#endif
