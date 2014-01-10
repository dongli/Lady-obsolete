#ifndef __Lady_AdvectionManager__
#define __Lady_AdvectionManager__

#include "TracerManager.h"

namespace lady {

class AdvectionManager {
protected:
    TracerManager tracerManager;
    LADY_REGRID *regrid; // used to interpolate velocity
public:
    AdvectionManager();
    ~AdvectionManager();

    /**
     *  Initialize advection manager.
     *
     *  @param domain    the spatial domain.
     *  @param mesh      the mesh.
     *  @param numParcel the number of parcels (a.k.a., tracers).
     */
    void init(const LADY_DOMAIN &domain, const LADY_MESH &mesh, int numParcel);

    /**
     *  Input one tracer species from given scalar field.
     *
     *  @param longName the long name of the tracer.
     *  @param q        the input scalar field.
     */
    void input(const string &longName, const LADY_SCALAR_FIELD &q);

    /**
     *  Output tracers on old time level into netCDF file.
     *
     *  @param fileName   the file name.
     *  @param newTimeIdx the new time level index.
     */
    void output(const string &fileName, const TimeLevelIndex<2> &newTimeIdx);

    /**
     *  Advect tracers one time step forward. 4th-order Runge-Kutta method is
     *  used to integrate advection equations.
     *
     *  @param dt         the time step size.
     *  @param oldTimeIdx the old time level index.
     *  @param V          the velocity field.
     */
    void advance(double dt, const TimeLevelIndex<2> &oldTimeIdx,
                 const LADY_VELOCITY_FIELD &V);};

}

#endif