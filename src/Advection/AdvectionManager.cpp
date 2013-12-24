#include "AdvectionManager.h"

namespace lady {

AdvectionManager::AdvectionManager() {
    regrid = NULL;
    REPORT_ONLINE;
}

AdvectionManager::~AdvectionManager() {
    if (regrid != NULL) {
        delete regrid;
    }
    REPORT_OFFLINE;
}

void AdvectionManager::init(const LADY_DOMAIN &domain, const LADY_MESH &mesh,
                            int numParcel) {
    tracerManager.init(domain, mesh, numParcel);
}

void AdvectionManager::output(const string &fileName, int timeLevel) {
    tracerManager.output(fileName, timeLevel);
}
    
/*
    Trajectory equation:

                                        dx
                                        -- = v,
                                        dt

    Deformation equation:

                                      dH
                                      -- = ∇v H.
                                      dt
 */

void AdvectionManager::advance(double dt, int oldLevel, int halfLevel,
                               int newLevel,
                               const LADY_VELOCITY_FIELD &V,
                               const LADY_TENSOR_FIELD &T) {
    const LADY_MESH &mesh = static_cast<const LADY_MESH&>(V.getMesh());
    const LADY_DOMAIN &domain = static_cast<const LADY_DOMAIN&>(mesh.getDomain());
    if (regrid == NULL) {
        regrid = new LADY_REGRID(mesh);
    }
    double dt05 = 0.5*dt;
    LADY_LIST<Tracer*>::iterator tracer = tracerManager.tracers.begin();
    for (; tracer != tracerManager.tracers.end(); ++tracer) {
        LADY_SPACE_COORD &x0 = (*tracer)->getX(oldLevel);
        LADY_SPACE_COORD &x1 = (*tracer)->getX(newLevel);
        LADY_MATRIX &H0 = (*tracer)->getH(oldLevel);
        LADY_MATRIX &H1 = (*tracer)->getH(newLevel);
        LADY_MESH_INDEX &idx0 = (*tracer)->getMeshIndex(oldLevel);
        LADY_MESH_INDEX &idx1 = (*tracer)->getMeshIndex(newLevel);
        LADY_VELOCITY v1(domain.getNumDim());
        LADY_VELOCITY v2(domain.getNumDim());
        LADY_VELOCITY v3(domain.getNumDim());
        LADY_VELOCITY v4(domain.getNumDim());
        LADY_VELOCITY v(domain.getNumDim());
        LADY_MATRIX t(domain.getNumDim(), domain.getNumDim());
        LADY_MATRIX dH1(domain.getNumDim(), domain.getNumDim());
        LADY_MATRIX dH2(domain.getNumDim(), domain.getNumDim());
        LADY_MATRIX dH3(domain.getNumDim(), domain.getNumDim());
        LADY_MATRIX dH4(domain.getNumDim(), domain.getNumDim());
        // TODO: Should we hide the following codes? Because they are related to
        //       sphere domain.
        if (idx0.isOnPole()) {
            idx0.setMoveOnPole(true);
            idx1.setMoveOnPole(true);
        } else {
            idx0.setMoveOnPole(false);
            idx1.setMoveOnPole(false);
        }
        regrid->run(BILINEAR, oldLevel, V, x0, v1, &idx0);
        regrid->run(BILINEAR, oldLevel, T, x0, t, &idx0); dH1 = t*H0;
        // ---------------------------------------------------------------------
        // stage 1
        mesh.move(x0, dt05, v1, idx0, x1); idx1.locate(mesh, x1);
        H1 = H0+dt05*dH1;
        regrid->run(BILINEAR, halfLevel, V, x1, v2, &idx1);
        regrid->run(BILINEAR, halfLevel, T, x1, t, &idx1); dH2 = t*H1;
        // ---------------------------------------------------------------------
        // stage 2
        mesh.move(x0, dt05, v2, idx0, x1); idx1.locate(mesh, x1);
        H1 = H0+dt05*dH2;
        regrid->run(BILINEAR, halfLevel, V, x1, v3, &idx1);
        regrid->run(BILINEAR, halfLevel, T, x1, t, &idx1); dH3 = t*H1;
        // ---------------------------------------------------------------------
        // stage 3
        mesh.move(x0, dt, v3, idx0, x1); idx1.locate(mesh, x1);
        H1 = H0+dt*dH3;
        regrid->run(BILINEAR, newLevel, V, x1, v4, &idx1);
        regrid->run(BILINEAR, newLevel, T, x1, t, &idx1); dH4 = t*H1;
        // ---------------------------------------------------------------------
        // stage 4
        v = (v1+v2*2.0+v3*2.0+v4)/6.0;
        mesh.move(x0, dt, v, idx0, x1); idx1.locate(mesh, x1);
        H1 = H0+dt*(dH1+dH2*2.0+dH3*2.0+dH4)/6.0;
    }
}

}
