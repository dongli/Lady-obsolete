#include "AdvectionManager.h"
#include "ShapeFunction.h"
#include "TracerSkeleton.h"

namespace lady {

AdvectionManager::AdvectionManager() {
    regrid = NULL;
    tracerMeshCells = NULL;
    REPORT_ONLINE;
}

AdvectionManager::~AdvectionManager() {
    if (regrid != NULL) {
        delete regrid;
    }
    if (tracerMeshCells != NULL) {
        delete tracerMeshCells;
    }
    REPORT_OFFLINE;
}

void AdvectionManager::init(const LADY_DOMAIN &domain, const LADY_MESH &mesh,
                            int numParcel) {
    tracerManager.init(domain, mesh, numParcel);
    ShapeFunction::init(domain);
    if (regrid == NULL) {
        regrid = new LADY_REGRID(mesh);
    }
    if (tracerMeshCells == NULL) {
        tracerMeshCells = new LADY_FIELD<TracerMeshCell>(mesh);
        tracerMeshCells->create(ScalarField, domain.getNumDim(), A_GRID);
        for (int i = 0; i < mesh.getTotalNumGrid(A_GRID); ++i) {
            LADY_SPACE_COORD x(domain.getNumDim());
            mesh.getGridCoord(i, x, A_GRID);
            x.transformToPS(domain);
            double volume;
            mesh.getCellVolume(i, volume);
            for (int l = 0; l < 2; ++l) {
                (*tracerMeshCells)(l, i).setCoord(x);
                (*tracerMeshCells)(l, i).setVolume(volume);
            }
        }
        cellCoords.reshape(3, mesh.getTotalNumGrid(A_GRID));
        for (int i = 0; i < mesh.getTotalNumGrid(A_GRID); ++i) {
            LADY_SPACE_COORD x = (*tracerMeshCells)(0, i).getCoord();
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            // TODO: How to hide the mesh details? The transformation should be
            //       moved into other low level place, and every coordinate
            //       class should provide getCartCoord() method.
            x.transformToCart(domain);
            cellCoords.col(i) = x.getCartCoord();
            // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        }
        cellTree = new Tree(cellCoords, cellCoordsMap);
        cellCoords = cellTree->Dataset();
    }
    TimeLevelIndex<2> timeIdx;
    assert(timeIdx.get() == 0);
    connectTracersAndMesh(timeIdx);
}

void AdvectionManager::registerTracer(const string &name, const string &units,
                                      const string &brief) {
    tracerManager.registerTracer(name, units, brief);
    for (int l = 0; l < 2; ++l) {
        for (int i = 0; i < tracerMeshCells->getMesh().getTotalNumGrid(A_GRID); ++i) {
            (*tracerMeshCells)(l, i).addSpecies();
        }
    }
}

void AdvectionManager::input(const string &name, const LADY_SCALAR_FIELD &q) {
    int speciesIdx = tracerManager.getSpeciesIndex(name);
    TimeLevelIndex<2> timeIdx;
    assert(timeIdx.get() == 0);
    // -------------------------------------------------------------------------
    // transfer the input tracer density field into internal tracer mass field
    assert(q.getMesh().getTotalNumGrid(A_GRID) ==
           tracerMeshCells->getMesh().getTotalNumGrid(A_GRID));
    for (int i = 0; i < q.getMesh().getTotalNumGrid(A_GRID); ++i) {
        double &m = (*tracerMeshCells)(timeIdx, i).getSpeciesMass(speciesIdx);
        m = q(timeIdx, i)*(*tracerMeshCells)(timeIdx, i).getVolume();
    }
    // -------------------------------------------------------------------------
    // transfer the tracer mass from cells to tracers
    LADY_LIST<Tracer*>::iterator tracer = tracerManager.tracers.begin();
    for (; tracer != tracerManager.tracers.end(); ++tracer) {
        double &m = (*tracer)->getSpeciesMass(timeIdx, speciesIdx);
        const list<TracerMeshCell*> &cells = (*tracer)->getConnectedCells();
        list<TracerMeshCell*>::const_iterator cell;
        for (cell = cells.begin(); cell != cells.end(); ++cell) {
            double weight = (*cell)->getWeight(*tracer)/(*cell)->getTotalRemapWeight();
            m += (*cell)->getSpeciesMass(speciesIdx)*weight;
        }
    }
    
}

void AdvectionManager::output(const string &fileName,
                              const TimeLevelIndex<2> &oldTimeIdx) {
    tracerManager.output(fileName, oldTimeIdx);
    // output the tracer density on the mesh
    
}

void AdvectionManager::advance(double dt, const TimeLevelIndex<2> &newTimeIdx,
                               const geomtk::RLLVelocityField &V) {
    integrate_RK4(dt, newTimeIdx, V);
    connectTracersAndMesh(newTimeIdx);
}

// -----------------------------------------------------------------------------
// private member functions

/*
    Trajectory equation:

                                        dx
                                        -- = v,
                                        dt

    Deformation equation (H is not computed from this equation):

                                      dH
                                      -- = ∇v H.
                                      dt
 */

void AdvectionManager::integrate_RK4(double dt,
                                     const TimeLevelIndex<2> &newTimeIdx,
                                     const LADY_VELOCITY_FIELD &V) {
    TimeLevelIndex<2> oldTimeIdx = newTimeIdx-1;
    TimeLevelIndex<2> halfTimeIdx = newTimeIdx-0.5;
    const LADY_MESH &mesh = static_cast<const LADY_MESH&>(V.getMesh());
    const LADY_DOMAIN &domain = static_cast<const LADY_DOMAIN&>(mesh.getDomain());
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

void AdvectionManager::connectTracersAndMesh(const TimeLevelIndex<2> &timeIdx) {
    const LADY_MESH &mesh = static_cast<const LADY_MESH&>(tracerMeshCells->getMesh());
    const LADY_DOMAIN &domain = static_cast<const LADY_DOMAIN&>(mesh.getDomain());
    // -------------------------------------------------------------------------
    // reset the connection between tracers and cells
    LADY_LIST<Tracer*>::iterator tracer = tracerManager.tracers.begin();
    for (; tracer != tracerManager.tracers.end(); ++tracer) {
        (*tracer)->resetConnectedCells();
    }
    for (int i = 0; i < mesh.getTotalNumGrid(A_GRID); ++i) {
        (*tracerMeshCells)(timeIdx, i).resetConnectedTracers();
    }
    // -------------------------------------------------------------------------
    // call mlpack::range::RangeSearch to find out the neighbor cells of tracers
    // and set the data structures for both cells and tracers for remapping
    tracer = tracerManager.tracers.begin();
    for (; tracer != tracerManager.tracers.end(); ++tracer) {
        // =====================================================================
        // search neighbor cells for the tracer
        LADY_SPACE_COORD &x = (*tracer)->getX(timeIdx);
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        // TODO: How to hide the mesh details? The transformation should be
        //       moved into other low level place, and every coordinate
        //       class should provide getCartCoord() method.
        x.transformToCart(domain);
        Searcher a(cellTree, NULL, cellCoords, x.getCartCoord(), true);
        // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        double longAxisSize = (*tracer)->getShapeSize(timeIdx)(0);
        mlpack::math::Range r(0.0, longAxisSize);
        vector<vector<size_t> > neighbors;
        vector<vector<double> > distances;
        a.Search(r, neighbors, distances);
        // =====================================================================
        // set data structures for the tracer and its neighbor cells
        for (int i = 0; i < neighbors[0].size(); ++i) {
            int cellIdx = cellCoordsMap[neighbors[0][i]];
            TracerMeshCell *cell = &(*tracerMeshCells)(timeIdx, cellIdx);
            // calculate the tracer shape function for the cell
            LADY_BODY_COORD y(domain.getNumDim());
            (*tracer)->getBodyCoord(domain, timeIdx, cell->getCoord(), y);
            double f = (*tracer)->getShapeFunction(timeIdx, y);
            if (f > 0.0) {
                cell->connect(*tracer, f);
                (*tracer)->connect(cell);
            }
        }
#define CHECK_NEIGHBORS 0
#if CHECK_NEIGHBORS == 1
        std::ofstream file;
        file.open("neighbors.txt");
        file << "p0 = (/" << (*tracer)->getX(timeIdx)(0) << "," << (*tracer)->getX(timeIdx)(1) << "/)" << endl;
        file << "ngb = new((/" << neighbors[0].size() << ",2/), double)" << endl;
        for (int m = 0; m < 2; ++m) {
            file << "ngb(:," << m << ") = (/";
            for (int i = 0; i < neighbors.size(); ++i) {
                for (int j = 0; j < neighbors[i].size()-1; ++j) {
                    file << (*tracerMeshCells)(0, cellCoordsMap[neighbors[i][j]]).getCoord()(m) << ",";
                }
                file << (*tracerMeshCells)(0, cellCoordsMap[neighbors[i][neighbors[i].size()-1]]).getCoord()(m) << "/)" << endl;
            }
        }
        int n = 100;
        file << "c0 = new((/" << n << ",2/), double)" << endl;
        double dtheta = PI2/(n-1);
        double lat = M_PI_2-longAxisSize/domain.getRadius();
        vector<LADY_SPACE_COORD*> c0(n);
        for (int i = 0; i < n; ++i) {
            c0[i] = new LADY_SPACE_COORD(2);
            LADY_SPACE_COORD xr(2);
            double theta = i*dtheta;
            xr(0) = theta;
            xr(1) = lat;
            domain.rotateBack((*tracer)->getX(timeIdx), *(c0[i]), xr);
        }
        for (int m = 0; m < 2; ++m) {
            file << "c0(:," << m << ") = (/";
            for (int i = 0; i < c0.size()-1; ++i) {
                file << (*(c0[i]))(m) << ",";
            }
            file << (*(c0[c0.size()-1]))(m) << "/)" << endl;
        }
        for (int i = 0; i < n; ++i) {
            delete c0[i];
        }
        file.close();
        CHECK_POINT
#endif
    }
}

}
