#ifndef __Lady_TracerSkeleton__
#define __Lady_TracerSkeleton__

#include "lady_commons.h"
#include "Tracer.h"

namespace lady {

class TracerSkeleton {
protected:
    Tracer *host;
    TimeLevels<vector<LADY_SPACE_COORD*>, 2> x;
    vector<LADY_BODY_COORD*> y; // fixed body coordinates
    TimeLevels<vector<LADY_MESH_INDEX*>, 2> idx;
public:
    TracerSkeleton(Tracer *host, int numDim);
    virtual ~TracerSkeleton();

    TracerSkeleton& operator=(const TracerSkeleton &other);

    void init(const LADY_DOMAIN &domain, const LADY_MESH &mesh);

    vector<LADY_SPACE_COORD*>& getXs(const TimeLevelIndex<2> &timeIdx);
    vector<LADY_BODY_COORD*>& getYs();
    vector<LADY_MESH_INDEX*>& getIdxs(const TimeLevelIndex<2> &timeIdx);
};

}

#endif