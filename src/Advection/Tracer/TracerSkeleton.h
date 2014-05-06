#ifndef __Lady_TracerSkeleton__
#define __Lady_TracerSkeleton__

#include "lady_commons.h"
#include "Tracer.h"

namespace lady {

class TracerSkeleton {
protected:
    Tracer *host;
    TimeLevels<vector<LADY_SPACE_COORD*>, 2> x; //>! spatial coordinates
    vector<LADY_BODY_COORD*> y;                 //>! fixed body coordinates
    TimeLevels<vector<LADY_MESH_INDEX*>, 2> idx;
    // Note: In sphere domain, xl are local stereographic projection coordinates,
    //       and in normal Cartesian domain, they should be the same with x.
    TimeLevels<vector<vec>, 2> xl;              //>! local coordinates
public:
    TracerSkeleton(Tracer *host, int numDim);
    virtual ~TracerSkeleton();

    TracerSkeleton& operator=(const TracerSkeleton &other);

    void init(const LADY_DOMAIN &domain, const LADY_MESH &mesh, double size);

    void updateLocalCoord(const LADY_DOMAIN &domain,
                          const TimeLevelIndex<2> &timeIdx);

    vector<LADY_SPACE_COORD*>& getSpaceCoords(const TimeLevelIndex<2> &timeIdx) {
        return x.getLevel(timeIdx);
    }

    const vector<LADY_BODY_COORD*>& getBodyCoords() const { return y; }

    vector<LADY_MESH_INDEX*>& getMeshIdxs(const TimeLevelIndex<2> &timeIdx) {
        return idx.getLevel(timeIdx);
    }

    const vector<vec>& getLocalCoords(const TimeLevelIndex<2> &timeIdx) const {
        return xl.getLevel(timeIdx);
    }
};

}

#endif
