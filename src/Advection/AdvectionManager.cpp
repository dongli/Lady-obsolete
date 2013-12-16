#include "AdvectionManager.h"

namespace lady {

AdvectionManager::AdvectionManager() {
    regrid = NULL;
}

AdvectionManager::~AdvectionManager() {
    if (regrid != NULL) {
        delete regrid;
    }
}

void AdvectionManager::init(const LADY_DOMAIN &domain, const LADY_MESH &mesh,
                            int numParcel) {
    tracerManager.init(domain, mesh, numParcel);
}

void AdvectionManager::output(const string &fileName, int timeLevel) {
    tracerManager.output(fileName, timeLevel);
}

void AdvectionManager::advance(double dt, int oldLevel, int halfLevel, int newLevel,
                               const LADY_VELOCITY_FIELD &V) {
    const LADY_MESH &mesh = static_cast<const LADY_MESH&>(V.getMesh());
    const LADY_DOMAIN &domain = static_cast<const LADY_DOMAIN&>(mesh.getDomain());
    if (regrid == NULL) {
        regrid = new LADY_REGRID(mesh);
    }
    double dt05 = 0.5*dt;
    LADY_LIST<Tracer*>::iterator t = tracerManager.tracers.begin();
    for (; t != tracerManager.tracers.end(); ++t) {
        LADY_SPACE_COORD &x0 = (*t)->getX(oldLevel);
        LADY_SPACE_COORD &x1 = (*t)->getX(newLevel);
//        LADY_MATRIX &H0 = (*t)->getH(oldLevel);
//        LADY_MATRIX &H1 = (*t)->getH(newLevel);
        LADY_MESH_INDEX &idx0 = (*t)->getMeshIndex(oldLevel);
        LADY_MESH_INDEX &idx1 = (*t)->getMeshIndex(newLevel);
        LADY_VELOCITY v1(domain.getNumDim());
        LADY_VELOCITY v2(domain.getNumDim());
        LADY_VELOCITY v3(domain.getNumDim());
        LADY_VELOCITY v4(domain.getNumDim());
        LADY_VELOCITY v(domain.getNumDim());
        // TODO: Should we hide the following codes? Because they are related to
        //       sphere domain.
        if (idx0.isOnPole()) {
            idx0.toggleMoveOnPole();
            idx1.toggleMoveOnPole();
        }
        // ---------------------------------------------------------------------
        regrid->run(BILINEAR, oldLevel, V, x0, v1, &idx0);
        mesh.move(x0, dt05, v1, idx0, x1);
        idx1.locate(mesh, x1);
        regrid->run(BILINEAR, halfLevel, V, x1, v2, &idx1);
        // ---------------------------------------------------------------------
        mesh.move(x0, dt05, v2, idx0, x1);
        idx1.locate(mesh, x1);
        regrid->run(BILINEAR, halfLevel, V, x1, v3, &idx1);
        // ---------------------------------------------------------------------
        mesh.move(x0, dt, v3, idx0, x1);
        idx1.locate(mesh, x1);
        regrid->run(BILINEAR, newLevel, V, x1, v4, &idx1);
        // ---------------------------------------------------------------------
        v = (v1+v2*2.0+v3*2.0+v4)/6.0;
        mesh.move(x0, dt, v, idx0, x1);
        idx1.locate(mesh, x1);
    }
}

}
