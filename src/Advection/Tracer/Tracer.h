#ifndef __Lady_Tracer__
#define __Lady_Tracer__

#include "commons.h"

namespace lady {

/**
 * Tracer
 * This class describes the tracer that is used to be advected by external wind
 * field.
 */

class Tracer {
protected:
    TimeLevels<LADY_SPACE_COORD*, 2> q; // centroid coordinate
    TimeLevels<LADY_MATRIX*, 2> H; // deformation matrix
public:
    Tracer(int numDim);
    virtual ~Tracer();

    Tracer& operator=(const Tracer &other);

    LADY_SPACE_COORD& getX(int timeLevel);
    LADY_MATRIX& getH(int timeLevel);
};

}

#endif
