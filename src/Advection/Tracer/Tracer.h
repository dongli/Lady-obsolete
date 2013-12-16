#ifndef __Lady_Tracer__
#define __Lady_Tracer__

#include "lady_commons.h"

namespace lady {

/**
 * Tracer
 * This class describes the tracer that is used to be advected by external wind
 * field.
 */

class Tracer {
protected:
    int ID;
    geomtk::TimeLevels<LADY_SPACE_COORD*, 2> q; // centroid coordinate
    geomtk::TimeLevels<LADY_MATRIX*, 2> H; // deformation matrix
    geomtk::TimeLevels<LADY_MESH_INDEX*, 2> idx; // tracer mesh index
public:
    Tracer(int numDim);
    virtual ~Tracer();

    Tracer& operator=(const Tracer &other);

    int getID() const;
    void setID(int ID);
    LADY_SPACE_COORD& getX(int timeLevel);
    LADY_MATRIX& getH(int timeLevel);
    LADY_MESH_INDEX& getMeshIndex(int timeLevel);
};

}

#endif
