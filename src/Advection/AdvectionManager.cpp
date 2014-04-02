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
    delete cellTree;
    REPORT_OFFLINE;
}

void AdvectionManager::init(const LADY_DOMAIN &domain, const LADY_MESH &mesh,
                            int numParcelX, int numParcelY) {
    TimeLevelIndex<2> initTimeIdx;
    assert(initTimeIdx.get() == 0);
    // -------------------------------------------------------------------------
    // initialize tracer manager
    tracerManager.init(domain, mesh, numParcelX, numParcelY);
    // -------------------------------------------------------------------------
    // initialize shape function (or kernel function)
    ShapeFunction::init(domain);
    // -------------------------------------------------------------------------
    // initialize regrid object
    if (regrid == NULL) {
        regrid = new LADY_REGRID(mesh);
    }
    // -------------------------------------------------------------------------
    // initialize tracer mesh cells
    tracerMeshCells.create("", "", "mesh cells for storing tracers", mesh, CENTER);
    for (int i = 0; i < mesh.getTotalNumGrid(CENTER); ++i) {
        LADY_SPACE_COORD x(domain.getNumDim());
        mesh.getGridCoord(i, CENTER, x);
        x.transformToCart(domain);
        double volume = mesh.getCellVolume(i);
        for (int l = 0; l < 2; ++l) {
            tracerMeshCells(initTimeIdx+l, i).setCoord(x);
            tracerMeshCells(initTimeIdx+l, i).setVolume(volume);
            tracerMeshCells(initTimeIdx+l, i).setID(i);
        }
    }
    // -------------------------------------------------------------------------
    // initialize tree structure of mesh grids
    cellCoords.reshape(3, mesh.getTotalNumGrid(CENTER));
    for (int i = 0; i < mesh.getTotalNumGrid(CENTER); ++i) {
        cellCoords.col(i) = tracerMeshCells(initTimeIdx, i).getCoord().getCartCoord();
    }
    cellTree = new Tree(cellCoords, cellCoordsMap);
    cellCoords = cellTree->Dataset();
    // -------------------------------------------------------------------------
    connectTracersAndMesh(initTimeIdx);
}

void AdvectionManager::registerTracer(const string &name, const string &units,
                                      const string &brief) {
    tracerManager.registerTracer(name, units, brief);
    TimeLevelIndex<2> initTimeIdx;
    for (int l = 0; l < 2; ++l) {
        totalMass.getLevel(l).push_back(0);
        for (int i = 0; i < tracerMeshCells.getMesh().getTotalNumGrid(CENTER); ++i) {
            tracerMeshCells(initTimeIdx+l, i).addSpecies();
        }
    }
}

void AdvectionManager::input(const TimeLevelIndex<2> &timeIdx,
                             vector<LADY_SCALAR_FIELD*> &q) {
    assert(timeIdx.isCurrentIndex());
    assert(q.size() == tracerManager.getNumSpecies());
    const LADY_MESH &mesh = static_cast<const LADY_MESH&>(tracerMeshCells.getMesh());
    // -------------------------------------------------------------------------
    // copy the input tracer density onto internal mesh grids
    for (int s = 0; s < q.size(); ++s) {
#ifdef DEBUG
        assert(q[s]->getMesh().getTotalNumGrid(CENTER) ==
               tracerMeshCells.getMesh().getTotalNumGrid(CENTER));
#endif
        for (int i = 0; i < mesh.getTotalNumGrid(CENTER); ++i) {
            double &rho = tracerMeshCells(timeIdx, i).getSpeciesDensity(s);
            rho = (*q[s])(timeIdx, i);
        }
        calcTotalMass(timeIdx);
    }
    // -------------------------------------------------------------------------
    // transfer the tracer mass from cells to tracers
    remapMeshToTracers(timeIdx);
    diagnose(timeIdx);
    // TODO: The mass on cells could be different after remapping from tracers.
    remapTracersToMesh(timeIdx);
    diagnose(timeIdx);
}

void AdvectionManager::output(const string &fileName,
                              const TimeLevelIndex<2> &oldTimeIdx) {
    tracerManager.output(fileName, oldTimeIdx);
    // output the tracer density on the mesh
    const LADY_MESH &mesh = static_cast<const LADY_MESH&>(tracerMeshCells.getMesh());
    const LADY_DOMAIN &domain = static_cast<const LADY_DOMAIN&>(mesh.getDomain());
    int ncId, lonDimId, latDimId;
    int lonVarId, latVarId;
    int dimIds[domain.getNumDim()];
    int qVarIds[tracerManager.getNumSpecies()], volVarId;
    vec lon, lat;
    
    if (nc_open(fileName.c_str(), NC_WRITE, &ncId) != NC_NOERR) {
        REPORT_ERROR("Failed to open " << fileName << "!");
    }
    nc_redef(ncId);
    // -------------------------------------------------------------------------
    // define dimensions
    // =========================================================================
    // longitude dimension
    if (nc_def_dim(ncId, "lon", mesh.getNumGrid(0, FULL), &lonDimId)
        != NC_NOERR) {
        REPORT_ERROR("Failed to define dimension lon!");
    }
    if (nc_def_var(ncId, "lon", NC_DOUBLE, 1, &lonDimId, &lonVarId)
        != NC_NOERR) {
        REPORT_ERROR("Failed to define coordinate variable lon!");
    }
    nc_put_att(ncId, lonVarId, "long_name", NC_CHAR, 9, "longitude");
    nc_put_att(ncId, lonVarId, "units", NC_CHAR, 12, "degrees_east");
    // =========================================================================
    // latitude dimension
    if (nc_def_dim(ncId, "lat", mesh.getNumGrid(1, FULL), &latDimId)
        != NC_NOERR) {
        REPORT_ERROR("Failed to define dimension lat!");
    }
    if (nc_def_var(ncId, "lat", NC_DOUBLE, 1, &latDimId, &latVarId)
        != NC_NOERR) {
        REPORT_ERROR("Failed to define coordinate variable lat!");
    }
    nc_put_att(ncId, latVarId, "long_name", NC_CHAR, 8, "latitude");
    nc_put_att(ncId, latVarId, "units", NC_CHAR, 13, "degrees_north");
    // =========================================================================
    dimIds[0] = latDimId;
    dimIds[1] = lonDimId;
    for (int i = 0; i < tracerManager.getNumSpecies(); ++i) {
        const TracerSpeciesInfo &speciesInfo = tracerManager.getSpeciesInfo(i);
        if (nc_def_var(ncId, speciesInfo.getName().c_str(), NC_DOUBLE,
                       domain.getNumDim(), dimIds, &qVarIds[i])
            != NC_NOERR) {
            REPORT_ERROR("Failed to define variable " <<
                         speciesInfo.getName().c_str() << "!");
        }
        nc_put_att(ncId, qVarIds[i], "long_name", NC_CHAR,
                   speciesInfo.getBrief().length(), speciesInfo.getBrief().c_str());
        nc_put_att(ncId, qVarIds[i], "units", NC_CHAR,
                   speciesInfo.getUnits().length(), speciesInfo.getUnits().c_str());
    }
    if (nc_def_var(ncId, "volume", NC_DOUBLE, domain.getNumDim(), dimIds, &volVarId)
        != NC_NOERR) {
        REPORT_ERROR("Failed to define variable volume!");
    }
    nc_put_att(ncId, volVarId, "long_name", NC_CHAR, 16, "mesh cell volume");
    nc_enddef(ncId);
    // -------------------------------------------------------------------------
    // put variables
    // =========================================================================
    lon = mesh.getGridCoords(0, FULL)/RAD;
    lat = mesh.getGridCoords(1, FULL)/RAD;
    nc_put_var(ncId, lonVarId, lon.memptr());
    nc_put_var(ncId, latVarId, lat.memptr());
    double *x  = new double[mesh.getTotalNumGrid(CENTER)];
    for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
        for (int i = 0; i < mesh.getTotalNumGrid(CENTER); ++i) {
            x[i] = tracerMeshCells(oldTimeIdx, i).getSpeciesDensity(s);
        }
        nc_put_var(ncId, qVarIds[s], x);
    }
    for (int i = 0; i < mesh.getTotalNumGrid(CENTER); ++i) {
        x[i] = tracerMeshCells(oldTimeIdx, i).getVolume();
    }
    nc_put_var(ncId, volVarId, x);
    delete [] x;
    nc_close(ncId);
}

void AdvectionManager::diagnose(const TimeLevelIndex<2> &timeIdx) {
    calcTotalMass(timeIdx);
    // print total mass for each species
    for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
        REPORT_NOTICE("Total mass of \"" <<
                      tracerManager.getSpeciesInfo(s).getName() << "\" is " <<
                      setw(30) << setprecision(20) <<
                      totalMass.getLevel(timeIdx)[s] << ".");
    }
}

void AdvectionManager::advance(double dt, const TimeLevelIndex<2> &newTimeIdx,
                               const geomtk::RLLVelocityField &velocity) {
    static clock_t time1, time2;
    
    time1 = clock();
    integrate_RK4(dt, newTimeIdx, velocity);
    time2 = clock();
    REPORT_NOTICE("integrate_RK uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");
    
    time1 = clock();
    connectTracersAndMesh(newTimeIdx);
    time2 = clock();
    REPORT_NOTICE("connectTracersAndMesh uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");
    
    time1 = clock();
    remapTracersToMesh(newTimeIdx);
    time2 = clock();
    REPORT_NOTICE("remapTracersToMesh uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");

//    diagnose(newTimeIdx);
}

// -----------------------------------------------------------------------------
// private member functions

void AdvectionManager::calcTotalMass(const TimeLevelIndex<2> &timeIdx) {
    // print total mass for each species
    for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
        totalMass.getLevel(timeIdx)[s] = 0;
        for (int i = 0; i < tracerMeshCells.getMesh().getTotalNumGrid(CENTER); ++i) {
            totalMass.getLevel(timeIdx)[s] +=
                tracerMeshCells(timeIdx, i).getSpeciesDensity(s)*
                tracerMeshCells(timeIdx, i).getVolume();
        }
    }
}

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
                                     const LADY_VELOCITY_FIELD &velocity) {
    TimeLevelIndex<2> oldTimeIdx = newTimeIdx-1;
    TimeLevelIndex<2> halfTimeIdx = newTimeIdx-0.5;
    const LADY_MESH &mesh = static_cast<const LADY_MESH&>(velocity.getMesh());
    const LADY_DOMAIN &domain = static_cast<const LADY_DOMAIN&>(mesh.getDomain());
    double dt05 = 0.5*dt;
#pragma omp parallel
#pragma omp single
    {
        LADY_LIST<Tracer*>::iterator tracer = tracerManager.tracers.begin();
        for (; tracer != tracerManager.tracers.end(); ++tracer) {
#pragma omp task
            {
                LADY_VELOCITY v1(domain.getNumDim());
                LADY_VELOCITY v2(domain.getNumDim());
                LADY_VELOCITY v3(domain.getNumDim());
                LADY_VELOCITY v4(domain.getNumDim());
                LADY_VELOCITY v(domain.getNumDim());
                const LADY_VELOCITY_FIELD::FieldType &divergence = velocity.getDivergence();
                double div;
                double k1_rho[tracerManager.getNumSpecies()];
                double k2_rho[tracerManager.getNumSpecies()];
                double k3_rho[tracerManager.getNumSpecies()];
                double k4_rho[tracerManager.getNumSpecies()];
                // -------------------------------------------------------------
                // update location and deformation matrix of tracer
                LADY_SPACE_COORD &x0 = (*tracer)->getX(oldTimeIdx);
                LADY_SPACE_COORD &x1 = (*tracer)->getX(newTimeIdx);
                LADY_MESH_INDEX &idx0 = (*tracer)->getMeshIndex(oldTimeIdx);
                LADY_MESH_INDEX &idx1 = (*tracer)->getMeshIndex(newTimeIdx);
                // TODO: Should we hide the following codes? Because they are
                //       related to sphere domain.
                if (idx0.isOnPole()) {
                    idx0.setMoveOnPole(true);
                    idx1.setMoveOnPole(true);
                    x0.transformToPS(domain);
                } else {
                    idx0.setMoveOnPole(false);
                    idx1.setMoveOnPole(false);
                }
                // =============================================================
                // stage 1
                regrid->run(BILINEAR, oldTimeIdx, velocity, x0, v1, &idx0);
                regrid->run(BILINEAR, oldTimeIdx, divergence, x0, div, &idx0);
                for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                    k1_rho[s] = -(*tracer)->getSpeciesDensity(s)*div;
                }
                mesh.move(x0, dt05, v1, idx0, x1); idx1.locate(mesh, x1);
                // =============================================================
                // stage 2
                regrid->run(BILINEAR, halfTimeIdx, velocity, x1, v2, &idx1);
                regrid->run(BILINEAR, halfTimeIdx, divergence, x1, div, &idx1);
                for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                    k2_rho[s] = -(*tracer)->getSpeciesDensity(s)*div;
                }
                mesh.move(x0, dt05, v2, idx0, x1); idx1.locate(mesh, x1);
                // =============================================================
                // stage 3
                regrid->run(BILINEAR, halfTimeIdx, velocity, x1, v3, &idx1);
                regrid->run(BILINEAR, newTimeIdx, divergence, x1, div, &idx1);
                for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                    k3_rho[s] = -(*tracer)->getSpeciesDensity(s)*div;
                }
                mesh.move(x0, dt, v3, idx0, x1); idx1.locate(mesh, x1);
                // =============================================================
                // stage 4
                regrid->run(BILINEAR, newTimeIdx, velocity, x1, v4, &idx1);
                regrid->run(BILINEAR, newTimeIdx, divergence, x1, div, &idx1);
                for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                    k4_rho[s] = -(*tracer)->getSpeciesDensity(s)*div;
                    (*tracer)->getSpeciesDensity(s) += dt*
                        (k1_rho[s]+2*k2_rho[s]+2*k3_rho[s]+k4_rho[s])/6;
                }
                v = (v1+v2*2.0+v3*2.0+v4)/6.0;
                mesh.move(x0, dt, v, idx0, x1); idx1.locate(mesh, x1);
                x1.transformToCart(domain);
                // -------------------------------------------------------------
                // update skeleton points of tracer
                TracerSkeleton &s = (*tracer)->getSkeleton();
                vector<LADY_SPACE_COORD*> &x0s = s.getSpaceCoords(oldTimeIdx);
                vector<LADY_SPACE_COORD*> &x1s = s.getSpaceCoords(newTimeIdx);
                vector<LADY_MESH_INDEX*> &idx0s = s.getMeshIdxs(oldTimeIdx);
                vector<LADY_MESH_INDEX*> &idx1s = s.getMeshIdxs(newTimeIdx);
                for (int i = 0; i < x0s.size(); ++i) {
                    if (idx0s[i]->isOnPole()) {
                        idx0s[i]->setMoveOnPole(true);
                        idx1s[i]->setMoveOnPole(true);
                        x0s[i]->transformToPS(domain);
                    } else {
                        idx0s[i]->setMoveOnPole(false);
                        idx1s[i]->setMoveOnPole(false);
                    }
                    // =========================================================
                    // stage 1
                    regrid->run(BILINEAR, oldTimeIdx, velocity, *x0s[i], v1, idx0s[i]);
                    mesh.move(*x0s[i], dt05, v1, *idx0s[i], *x1s[i]);
                    idx1s[i]->locate(mesh, *x1s[i]);
                    // =========================================================
                    // stage 2
                    regrid->run(BILINEAR, halfTimeIdx, velocity, *x1s[i], v2, idx1s[i]);
                    mesh.move(*x0s[i], dt05, v2, *idx0s[i], *x1s[i]);
                    idx1s[i]->locate(mesh, *x1s[i]);
                    // =========================================================
                    // stage 3
                    regrid->run(BILINEAR, halfTimeIdx, velocity, *x1s[i], v3, idx1s[i]);
                    mesh.move(*x0s[i], dt, v3, *idx0s[i], *x1s[i]);
                    idx1s[i]->locate(mesh, *x1s[i]);
                    // =========================================================
                    // stage 4
                    regrid->run(BILINEAR, newTimeIdx, velocity, *x1s[i], v4, idx1s[i]);
                    v = (v1+v2*2.0+v3*2.0+v4)/6.0;
                    mesh.move(*x0s[i], dt, v, *idx0s[i], *x1s[i]);
                    idx1s[i]->locate(mesh, *x1s[i]);
                    x1s[i]->transformToCart(domain);
                }
                // -------------------------------------------------------------
                (*tracer)->updateDeformMatrix(domain, mesh, newTimeIdx);
            }
        }
    }
}

void AdvectionManager::connectTracersAndMesh(const TimeLevelIndex<2> &timeIdx) {
    const LADY_MESH &mesh = static_cast<const LADY_MESH&>(tracerMeshCells.getMesh());
    const LADY_DOMAIN &domain = static_cast<const LADY_DOMAIN&>(mesh.getDomain());
    // -------------------------------------------------------------------------
    // reset the connection between tracers and cells
    LADY_LIST<Tracer*>::iterator tracer = tracerManager.tracers.begin();
    for (; tracer != tracerManager.tracers.end(); ++tracer) {
        (*tracer)->resetConnectedCells();
    }
    for (int i = 0; i < mesh.getTotalNumGrid(CENTER); ++i) {
        tracerMeshCells(timeIdx, i).resetConnectedTracers();
    }
    // -------------------------------------------------------------------------
    // call mlpack::range::RangeSearch to find out the neighbor cells of tracers
    // and set the data structures for both cells and tracers for remapping
    clock_t time1 = 0, time2 = 0;
    tracer = tracerManager.tracers.begin();
    for (; tracer != tracerManager.tracers.end(); ++tracer) {
        // =====================================================================
        // search neighbor cells for the tracer
        LADY_SPACE_COORD &x = (*tracer)->getX(timeIdx);
        clock_t clock1 = clock();
        Searcher a(cellTree, NULL, cellCoords, x.getCartCoord(), true);
        double longAxisSize = (*tracer)->getShapeSize(timeIdx)(0);
        mlpack::math::Range r(0.0, longAxisSize);
        vector<vector<size_t> > neighbors;
        vector<vector<double> > distances;
        a.Search(r, neighbors, distances);
        clock_t clock2 = clock();
        time1 += clock2-clock1;
        // =====================================================================
        // set data structures for the tracer and its neighbor cells
        clock1 = clock();
        LADY_BODY_COORD y(domain.getNumDim());
        for (int i = 0; i < neighbors[0].size(); ++i) {
            int cellIdx = cellCoordsMap[neighbors[0][i]];
            TracerMeshCell *cell = &tracerMeshCells(timeIdx, cellIdx);
            // calculate the tracer shape function for the cell
            (*tracer)->getBodyCoord(domain, timeIdx, cell->getCoord(), y);
            double f = (*tracer)->getShapeFunction(timeIdx, y);
            if (f > 0.0) {
                cell->connect(*tracer, f);
                (*tracer)->connect(cell, f);
            }
        }
        clock2 = clock();
        time2 = clock2-clock1;
        assert((*tracer)->getConnectedCells().size() != 0);
#define CHECK_NEIGHBORS 0
#if CHECK_NEIGHBORS == 1
        std::ofstream file;
        file.open("neighbors.txt");
        file << "p0 = (/" << (*tracer)->getX(timeIdx)(0) << "," << (*tracer)->getX(timeIdx)(1) << "/)" << endl;
        file << "ngb = new((/" << (*tracer)->getConnectedCells().size() << ",2/), double)" << endl;
        for (int m = 0; m < 2; ++m) {
            file << "ngb(:," << m << ") = (/";
            vector<TracerMeshCell*>::iterator cell = (*tracer)->getConnectedCells().begin();
            int i = 0;
            for (; cell != (*tracer)->getConnectedCells().end(); ++cell) {
                i++;
                if (i != (*tracer)->getConnectedCells().size()) {
                    file << (*cell)->getCoord()(m) << ",";
                } else {
                    file << (*cell)->getCoord()(m) << "/)" << endl;
                }
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
            xr.setCoord(theta, lat);
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
    REPORT_NOTICE("Range search uses " << setw(6) << setprecision(2) << (double)(time1)/CLOCKS_PER_SEC << " seconds.");
    REPORT_NOTICE("Shape function uses " << setw(6) << setprecision(2) << (double)(time2)/CLOCKS_PER_SEC << " seconds.");
}

void AdvectionManager::remapMeshToTracers(const TimeLevelIndex<2> &timeIdx) {
    LADY_LIST<Tracer*>::iterator tracer = tracerManager.tracers.begin();
    for (; tracer != tracerManager.tracers.end(); ++tracer) {
        (*tracer)->resetSpecies();
        const vector<TracerMeshCell*> &cells = (*tracer)->getConnectedCells();
        assert(cells.size() != 0);
        for (int i = 0; i < (*tracer)->getNumConnectedCell(); ++i) {
            double weight = cells[i]->getRemapWeight(*tracer)/
                            (*tracer)->getTotalRemapWeight();
            for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                double &rho = (*tracer)->getSpeciesDensity(s);
                rho += cells[i]->getSpeciesDensity(s)*weight;
            }
        }
    }
}

void AdvectionManager::remapTracersToMesh(const TimeLevelIndex<2> &timeIdx) {
    for (int i = 0; i < tracerMeshCells.getMesh().getTotalNumGrid(CENTER); ++i) {
        tracerMeshCells(timeIdx, i).resetSpecies();
    }
    LADY_LIST<Tracer*>::iterator tracer = tracerManager.tracers.begin();
    for (; tracer != tracerManager.tracers.end(); ++tracer) {
        vector<TracerMeshCell*> &cells = (*tracer)->getConnectedCells();
        for (int i = 0; i < (*tracer)->getNumConnectedCell(); ++i) {
            double weight = cells[i]->getRemapWeight(*tracer)/
                            cells[i]->getTotalRemapWeight();
            for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                double &rho = cells[i]->getSpeciesDensity(s);
                rho += (*tracer)->getSpeciesDensity(s)*weight;
            }
        }
    }
    correctTotalMassOnMesh(timeIdx);
}

void AdvectionManager::correctTotalMassOnMesh(const TimeLevelIndex<2> &timeIdx) {
    double expectedTotalMass[tracerManager.getNumSpecies()];
    if (timeIdx.isCurrentIndex()) {
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            expectedTotalMass[s] = totalMass.getLevel(timeIdx)[s];
        }
    } else {
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            expectedTotalMass[s] = totalMass.getLevel(timeIdx-1)[s];
        }
    }
    calcTotalMass(timeIdx);
    for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
        double fixer = expectedTotalMass[s]/totalMass.getLevel(timeIdx)[s];
        for (int i = 0; i < tracerMeshCells.getMesh().getTotalNumGrid(CENTER); ++i) {
            tracerMeshCells(timeIdx, i).getSpeciesDensity(s) *= fixer;
        }
    }
    calcTotalMass(timeIdx);
}

}
