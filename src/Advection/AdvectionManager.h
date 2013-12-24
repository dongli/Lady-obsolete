#ifndef __Lady_AdvectionManager__
#define __Lady_AdvectionManager__

#include "TracerManager.h"

namespace lady {

class AdvectionManager {
protected:
    TracerManager tracerManager;
    LADY_REGRID *regrid;
public:
    AdvectionManager();
    ~AdvectionManager();

    void init(const LADY_DOMAIN &domain, const LADY_MESH &mesh, int numParcel);

    void output(const string &fileName, int timeLevel);

    /**
     * Advance one time step.
     * 4th-order Runge-Kutta method is used to integrate advection equations.
     */
    void advance(double dt, int oldLevel, int halfLevel, int newLevel,
                 const LADY_VELOCITY_FIELD &V, const LADY_TENSOR_FIELD &T);
};

}

#endif