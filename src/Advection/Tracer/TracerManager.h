#ifndef __Lady_TracerManager__
#define __Lady_TracerManager__

#include "Tracer.h"

namespace lady {

class TracerManager {
protected:
    LADY_LIST<Tracer*> tracers;
public:
    TracerManager();
    virtual ~TracerManager();

    void init(const LADY_DOMAIN &domain, int numParcel);
};

}

#endif
