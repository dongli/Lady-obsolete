#ifndef __lady_commons__
#define __lady_commons__

#include "geomtk.h"

#include <iostream>
#include <iomanip>
#include <list>
#include <string>
#include <ctime>

namespace lady {
    
using std::cout;
using std::endl;
using std::setw;
using std::string;

using arma::vec;
using arma::cube;

using geomtk::RAD;
using geomtk::PI2;
using geomtk::BILINEAR;
using geomtk::Time;
using geomtk::A_GRID;
using geomtk::C_GRID;
using geomtk::CENTER;
using geomtk::EDGE;
    
#define LADY_SPACE_COORD geomtk::SphereCoord
#define LADY_VELOCITY geomtk::SphereVelocity
#define LADY_MATRIX arma::mat
#define LADY_LIST std::list

#define LADY_DOMAIN geomtk::SphereDomain
#define LADY_MESH geomtk::RLLMesh
#define LADY_MESH_INDEX geomtk::RLLMeshIndex
#define LADY_VELOCITY_FIELD geomtk::RLLVelocityField
#define LADY_REGRID geomtk::RLLRegrid

}

#endif
