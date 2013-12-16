#ifndef __Lady_TracerManager__
#define __Lady_TracerManager__

#include "Tracer.h"

namespace lady {

/**
 * TracerManager
 * This class stores the tracer objects and manages the input and output of them.
 */

class TracerManager {
    friend class AdvectionManager;
protected:
    const LADY_DOMAIN *domain;
    LADY_LIST<Tracer*> tracers;
public:
    TracerManager();
    virtual ~TracerManager();

    /**
     * Initializer of tracer manager.
     */
    void init(const LADY_DOMAIN &domain, const LADY_MESH &mesh, int numParcel);

    /**
     * Output tracers on given time level into netCDF file.
     */
    void output(const string &fileName, int timeLevel);
};

}

#endif
