#include "AdvectionManager.h"
#include "ShapeFunction.h"
#include "TracerSkeleton.h"

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
    ShapeFunction::init(domain);
}
    
void AdvectionManager::input(const string &longName,
                             const LADY_SCALAR_FIELD &q) {
//    const LADY_MESH &mesh = static_cast<const LADY_MESH&>(q.getMesh());
//    const LADY_DOMAIN &domain = static_cast<const LADY_DOMAIN&>(mesh.getDomain());
//    LADY_SPACE_COORD x(domain.getNumDim());
//    LADY_BODY_COORD y(domain.getNumDim());
//    tracerManager.registerTracer(longName);
//    LADY_LIST<Tracer*>::iterator tracer = tracerManager.tracers.begin();
//    for (; tracer != tracerManager.tracers.end(); ++tracer) {
//        LADY_MESH_INDEX &idx = (*tracer)->getMeshIndex(0);
//        // check which cells are covered by the tracer
//        mesh.getGridCoord(CENTER, idx, x);
//        (*tracer)->getBodyCoord(domain, 0, x, y);
//        double f = (*tracer)->getShapeFunction(0, y);
//    }
}

void AdvectionManager::output(const string &fileName,
                              const TimeLevelIndex<2> &oldTimeIdx) {
    tracerManager.output(fileName, oldTimeIdx);
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

void AdvectionManager::advance(double dt, const TimeLevelIndex<2> &newTimeIdx,
                               const LADY_VELOCITY_FIELD &V) {
    TimeLevelIndex<2> oldTimeIdx = newTimeIdx-1;
    TimeLevelIndex<2> halfTimeIdx = newTimeIdx-0.5;
    const LADY_MESH &mesh = static_cast<const LADY_MESH&>(V.getMesh());
    const LADY_DOMAIN &domain = static_cast<const LADY_DOMAIN&>(mesh.getDomain());
    if (regrid == NULL) {
        regrid = new LADY_REGRID(mesh);
    }
    double dt05 = 0.5*dt;
    LADY_LIST<Tracer*>::iterator tracer = tracerManager.tracers.begin();
    for (; tracer != tracerManager.tracers.end(); ++tracer) {
        // ---------------------------------------------------------------------
        // update location and deformation matrix of tracer
        LADY_SPACE_COORD &x0 = (*tracer)->getX(oldTimeIdx);
        LADY_SPACE_COORD &x1 = (*tracer)->getX(newTimeIdx);
        LADY_MESH_INDEX &idx0 = (*tracer)->getMeshIndex(oldTimeIdx);
        LADY_MESH_INDEX &idx1 = (*tracer)->getMeshIndex(newTimeIdx);
        LADY_VELOCITY v1(domain.getNumDim());
        LADY_VELOCITY v2(domain.getNumDim());
        LADY_VELOCITY v3(domain.getNumDim());
        LADY_VELOCITY v4(domain.getNumDim());
        LADY_VELOCITY v(domain.getNumDim());
        // TODO: Should we hide the following codes? Because they are related to
        //       sphere domain.
        if (idx0.isOnPole()) {
            idx0.setMoveOnPole(true);
            idx1.setMoveOnPole(true);
        } else {
            idx0.setMoveOnPole(false);
            idx1.setMoveOnPole(false);
        }
        regrid->run(BILINEAR, oldTimeIdx, V, x0, v1, &idx0);
        // =====================================================================
        // stage 1
        mesh.move(x0, dt05, v1, idx0, x1); idx1.locate(mesh, x1);
        regrid->run(BILINEAR, halfTimeIdx, V, x1, v2, &idx1);
        // =====================================================================
        // stage 2
        mesh.move(x0, dt05, v2, idx0, x1); idx1.locate(mesh, x1);
        regrid->run(BILINEAR, halfTimeIdx, V, x1, v3, &idx1);
        // =====================================================================
        // stage 3
        mesh.move(x0, dt, v3, idx0, x1); idx1.locate(mesh, x1);
        regrid->run(BILINEAR, newTimeIdx, V, x1, v4, &idx1);
        // =====================================================================
        // stage 4
        v = (v1+v2*2.0+v3*2.0+v4)/6.0;
        mesh.move(x0, dt, v, idx0, x1); idx1.locate(mesh, x1);
        // ---------------------------------------------------------------------
        // update skeleton points of tracer
        TracerSkeleton &s = (*tracer)->getSkeleton();
        vector<LADY_SPACE_COORD*> x0s = s.getXs(oldTimeIdx);
        vector<LADY_SPACE_COORD*> x1s = s.getXs(newTimeIdx);
        vector<LADY_MESH_INDEX*> idx0s = s.getIdxs(oldTimeIdx);
        vector<LADY_MESH_INDEX*> idx1s = s.getIdxs(newTimeIdx);
        for (int i = 0; i < x0s.size(); ++i) {
            LADY_VELOCITY v1(domain.getNumDim());
            LADY_VELOCITY v2(domain.getNumDim());
            LADY_VELOCITY v3(domain.getNumDim());
            LADY_VELOCITY v4(domain.getNumDim());
            LADY_VELOCITY v(domain.getNumDim());
            if (idx0s[i]->isOnPole()) {
                idx0s[i]->setMoveOnPole(true);
                idx1s[i]->setMoveOnPole(true);
            } else {
                idx0s[i]->setMoveOnPole(false);
                idx1s[i]->setMoveOnPole(false);
            }
            regrid->run(BILINEAR, oldTimeIdx, V, *x0s[i], v1, idx0s[i]);
            // =================================================================
            // stage 1
            mesh.move(*x0s[i], dt05, v1, idx0, *x1s[i]);
            idx1s[i]->locate(mesh, *x1s[i]);
            regrid->run(BILINEAR, halfTimeIdx, V, *x1s[i], v2, idx1s[i]);
            // =================================================================
            // stage 2
            mesh.move(*x0s[i], dt05, v2, idx0, *x1s[i]);
            idx1s[i]->locate(mesh, *x1s[i]);
            regrid->run(BILINEAR, halfTimeIdx, V, *x1s[i], v3, idx1s[i]);
            // =================================================================
            // stage 3
            mesh.move(*x0s[i], dt, v3, idx0, *x1s[i]);
            idx1s[i]->locate(mesh, *x1s[i]);
            regrid->run(BILINEAR, newTimeIdx, V, *x1s[i], v4, idx1s[i]);
            // =================================================================
            // stage 4
            v = (v1+v2*2.0+v3*2.0+v4)/6.0;
            mesh.move(*x0s[i], dt, v, *idx0s[i], *x1s[i]);
            idx1s[i]->locate(mesh, *x1s[i]);
        }
        // ---------------------------------------------------------------------
        (*tracer)->updateH(domain, newTimeIdx);
        (*tracer)->selfInspect(domain, newTimeIdx);
    }
}

}
