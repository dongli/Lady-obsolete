#ifndef __Lady_TracerManager__
#define __Lady_TracerManager__

#include "Tracer.h"
#include "TracerSpeciesInfo.h"

namespace lady {

/**
 *  This class stores the tracer objects and manages the initialization,
 *  registration and output.
 */
class TracerManager {
    friend class AdvectionManager;
protected:
    const LADY_DOMAIN *domain;
    LADY_LIST<Tracer*> tracers;
    vector<TracerSpeciesInfo*> speciesInfos;
public:
    TracerManager();
    virtual ~TracerManager();

    /**
     *  Initialize tracer manager.
     *
     *  @param domain    the space domain.
     *  @param mesh      the mesh where flow is defined.
     *  @param numParcel the number of parcels.
     */
    void init(const LADY_DOMAIN &domain, const LADY_MESH &mesh, int numParcel);

    void registerTracer(const string &name, const string &units,
                        const string &brief);

    int getSpeciesIndex(const string &name) const;

    /**
     *  Output tracers on old time level into netCDF file.
     *
     *  @param fileName   the output netCDF file name.
     *  @param oldTimeIdx the old time level index.
     */
    void output(const string &fileName, const TimeLevelIndex<2> &oldTimeIdx);
};

}

#endif
