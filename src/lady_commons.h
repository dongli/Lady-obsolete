#ifndef __lady_commons__
#define __lady_commons__

#include "geomtk.h"
#include "mlpack/methods/range_search/range_search.hpp"
#include "nlopt.hpp"

#include <iostream>
#include <iomanip>
#include <list>
#include <string>
#include <ctime>
#include <fstream>

namespace lady {
    
using std::cout;
using std::endl;
using std::setw;
using std::setprecision;
using std::string;
using std::vector;
using std::list;

using arma::vec;
using arma::mat;
using arma::cube;

using geomtk::RAD;
using geomtk::PI2;
using geomtk::BILINEAR;
using geomtk::Time;
using geomtk::ScalarField;
using geomtk::VectorField;
using geomtk::A_GRID;
using geomtk::C_GRID;
using geomtk::CENTER;
using geomtk::EDGE;
using geomtk::TimeLevels;
using geomtk::TimeLevelIndex;
using geomtk::SystemTools;
    
#define LADY_SPACE_COORD geomtk::SphereCoord
#define LADY_BODY_COORD geomtk::BodyCoord
#define LADY_VELOCITY geomtk::SphereVelocity
#define LADY_MATRIX mat
#define LADY_LIST list

#define LADY_DOMAIN geomtk::SphereDomain
#define LADY_MESH geomtk::RLLMesh
#define LADY_MESH_INDEX geomtk::RLLMeshIndex
#define LADY_FIELD geomtk::RLLField
#define LADY_SCALAR_FIELD geomtk::RLLField<double>
#define LADY_VELOCITY_FIELD geomtk::RLLVelocityField
#define LADY_REGRID geomtk::RLLRegrid

}

#endif
