#include "AdvectionManager.h"
#include "ShapeFunction.h"
#include "TracerSkeleton.h"

namespace lady {

AdvectionManager::AdvectionManager() {
    domain = NULL;
    mesh = NULL;
    regrid = NULL;
    numLongTracer = 0;
    numUnresolvedTracer = 0;
    numVoidCell = 0;
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
    this->domain = &domain;
    this->mesh = &mesh;
    TimeLevelIndex<2> initTimeIdx;
#ifndef NDEBUG
    assert(initTimeIdx.get() == 0);
#endif
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
        for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
            tracerMeshCells(initTimeIdx+l, i).addSpecies();
        }
    }
}

void AdvectionManager::input(const TimeLevelIndex<2> &timeIdx,
                             vector<LADY_SCALAR_FIELD*> &q) {
#ifndef NDEBUG
    assert(q.size() == tracerManager.getNumSpecies());
#endif
    // -------------------------------------------------------------------------
    // copy the input tracer density onto internal mesh grids
    for (int s = 0; s < q.size(); ++s) {
#ifndef NDEBUG
        assert(q[s]->getMesh().getTotalNumGrid(CENTER) ==
               tracerMeshCells.getMesh().getTotalNumGrid(CENTER));
#endif
        for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
            TracerMeshCell &cell = tracerMeshCells(timeIdx, i);
            cell.getSpeciesDensity(s) = (*q[s])(timeIdx, i);
            cell.getSpeciesMass(s) = (*q[s])(timeIdx, i)*cell.getVolume();
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
    int ncId, lonDimId, latDimId;
    int lonVarId, latVarId;
    int dimIds[domain->getNumDim()];
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
    if (nc_def_dim(ncId, "lon", mesh->getNumGrid(0, FULL), &lonDimId)
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
    if (nc_def_dim(ncId, "lat", mesh->getNumGrid(1, FULL), &latDimId)
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
                       domain->getNumDim(), dimIds, &qVarIds[i])
            != NC_NOERR) {
            REPORT_ERROR("Failed to define variable " <<
                         speciesInfo.getName().c_str() << "!");
        }
        nc_put_att(ncId, qVarIds[i], "long_name", NC_CHAR,
                   speciesInfo.getBrief().length(), speciesInfo.getBrief().c_str());
        nc_put_att(ncId, qVarIds[i], "units", NC_CHAR,
                   speciesInfo.getUnits().length(), speciesInfo.getUnits().c_str());
    }
    if (nc_def_var(ncId, "volume", NC_DOUBLE, domain->getNumDim(), dimIds, &volVarId)
        != NC_NOERR) {
        REPORT_ERROR("Failed to define variable volume!");
    }
    nc_put_att(ncId, volVarId, "long_name", NC_CHAR, 16, "mesh cell volume");
    nc_enddef(ncId);
    // -------------------------------------------------------------------------
    // put variables
    // =========================================================================
    lon = mesh->getGridCoords(0, FULL)/RAD;
    lat = mesh->getGridCoords(1, FULL)/RAD;
    nc_put_var(ncId, lonVarId, lon.memptr());
    nc_put_var(ncId, latVarId, lat.memptr());
    double *x  = new double[mesh->getTotalNumGrid(CENTER)];
    for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
        for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
            x[i] = tracerMeshCells(oldTimeIdx, i).getSpeciesDensity(s);
        }
        nc_put_var(ncId, qVarIds[s], x);
    }
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
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
    embedTracersIntoMesh(newTimeIdx);
    time2 = clock();
    REPORT_NOTICE("embedTracersIntoMesh uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");

    time1 = clock();
    connectTracersAndMesh(newTimeIdx);
    time2 = clock();
    REPORT_NOTICE("connectTracersAndMesh uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");

    time1 = clock();
    splitTracers(newTimeIdx);
    time2 = clock();
    REPORT_NOTICE("splitTracers uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");

    time1 = clock();
    mergeTracers(newTimeIdx);
    time2 = clock();
    REPORT_NOTICE("mergeTracers uses " << setw(6) << setprecision(2) << (double)(time2-time1)/CLOCKS_PER_SEC << " seconds.");

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
        for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
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
    double dt05 = 0.5*dt;
    numLongTracer = 0;
#pragma omp parallel
#pragma omp single
    {
        LADY_LIST<Tracer*>::iterator tracer = tracerManager.tracers.begin();
        for (; tracer != tracerManager.tracers.end(); ++tracer) {
#pragma omp task
            {
                LADY_VELOCITY v1(domain->getNumDim());
                LADY_VELOCITY v2(domain->getNumDim());
                LADY_VELOCITY v3(domain->getNumDim());
                LADY_VELOCITY v4(domain->getNumDim());
                LADY_VELOCITY v(domain->getNumDim());
                const LADY_VELOCITY_FIELD::FieldType &divergence = velocity.getDivergence();
                double div;
                vec rho(tracerManager.tracers.size());
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
                    x0.transformToPS(*domain);
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
                    rho[s] = (*tracer)->getSpeciesDensity(s)+dt05*k1_rho[s];
                }
                mesh->move(x0, dt05, v1, idx0, x1);
                idx1.locate(*mesh, x1);
                // =============================================================
                // stage 2
                regrid->run(BILINEAR, halfTimeIdx, velocity, x1, v2, &idx1);
                regrid->run(BILINEAR, halfTimeIdx, divergence, x1, div, &idx1);
                for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                    k2_rho[s] = -rho[s]*div;
                    rho[s] = (*tracer)->getSpeciesDensity(s)+dt05*k2_rho[s];
                }
                mesh->move(x0, dt05, v2, idx0, x1);
                idx1.locate(*mesh, x1);
                // =============================================================
                // stage 3
                regrid->run(BILINEAR, halfTimeIdx, velocity, x1, v3, &idx1);
                regrid->run(BILINEAR, newTimeIdx, divergence, x1, div, &idx1);
                for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                    k3_rho[s] = -rho[s]*div;
                    rho[s] = (*tracer)->getSpeciesDensity(s)+dt*k3_rho[s];
                }
                mesh->move(x0, dt, v3, idx0, x1);
                idx1.locate(*mesh, x1);
                // =============================================================
                // stage 4
                regrid->run(BILINEAR, newTimeIdx, velocity, x1, v4, &idx1);
                regrid->run(BILINEAR, newTimeIdx, divergence, x1, div, &idx1);
                for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                    k4_rho[s] = -rho[s]*div;
                    (*tracer)->getSpeciesDensity(s) += dt*
                        (k1_rho[s]+2.0*k2_rho[s]+2.0*k3_rho[s]+k4_rho[s])/6.0;
                }
                v = (v1+v2*2.0+v3*2.0+v4)/6.0;
                mesh->move(x0, dt, v, idx0, x1);
                idx1.locate(*mesh, x1);
                x1.transformToCart(*domain);
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
                        x0s[i]->transformToPS(*domain);
                    } else {
                        idx0s[i]->setMoveOnPole(false);
                        idx1s[i]->setMoveOnPole(false);
                    }
                    // =========================================================
                    // stage 1
                    regrid->run(BILINEAR, oldTimeIdx, velocity, *x0s[i], v1, idx0s[i]);
                    mesh->move(*x0s[i], dt05, v1, *idx0s[i], *x1s[i]);
                    idx1s[i]->locate(*mesh, *x1s[i]);
                    // =========================================================
                    // stage 2
                    regrid->run(BILINEAR, halfTimeIdx, velocity, *x1s[i], v2, idx1s[i]);
                    mesh->move(*x0s[i], dt05, v2, *idx0s[i], *x1s[i]);
                    idx1s[i]->locate(*mesh, *x1s[i]);
                    // =========================================================
                    // stage 3
                    regrid->run(BILINEAR, halfTimeIdx, velocity, *x1s[i], v3, idx1s[i]);
                    mesh->move(*x0s[i], dt, v3, *idx0s[i], *x1s[i]);
                    idx1s[i]->locate(*mesh, *x1s[i]);
                    // =========================================================
                    // stage 4
                    regrid->run(BILINEAR, newTimeIdx, velocity, *x1s[i], v4, idx1s[i]);
                    v = (v1+v2*2.0+v3*2.0+v4)/6.0;
                    mesh->move(*x0s[i], dt, v, *idx0s[i], *x1s[i]);
                    idx1s[i]->locate(*mesh, *x1s[i]);
                    x1s[i]->transformToCart(*domain);
                }
                // -------------------------------------------------------------
                (*tracer)->updateDeformMatrix(*domain, *mesh, newTimeIdx);
                if ((*tracer)->getBadType() == Tracer::EXTREME_FILAMENTATION) {
                    if (numLongTracer == longTracers.size()) {
                        longTracers.push_back(tracer);
                    } else {
                        longTracers[numLongTracer] = tracer;
                    }
                    numLongTracer++;
                }
            }
        }
    }
}

void AdvectionManager::embedTracersIntoMesh(const TimeLevelIndex<2> &timeIdx) {
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
        tracerMeshCells(timeIdx, i).resetContainedTracers();
    }
    LADY_LIST<Tracer*>::iterator tracer = tracerManager.tracers.begin();
    for (; tracer != tracerManager.tracers.end(); ++tracer) {
        int i = mesh->wrapIndex((*tracer)->getMeshIndex(timeIdx), CENTER);
        TracerMeshCell &cell = tracerMeshCells(timeIdx, i);
        cell.contain(*tracer);
    }
}

void AdvectionManager::connectTracerAndMesh(const TimeLevelIndex<2> &timeIdx,
                                            list<Tracer*>::iterator &tracer) {
    // call mlpack::range::RangeSearch to find out the neighbor cells of tracers
    // and set the data structures for both cells and tracers for remapping
    Searcher a(cellTree, NULL, cellCoords,
               (*tracer)->getX(timeIdx).getCartCoord(), true);
    double longAxisSize = (*tracer)->getShapeSize(timeIdx)[0];
    mlpack::math::Range r(0.0, longAxisSize);
    vector<vector<size_t> > neighbors;
    vector<vector<double> > distances;
    a.Search(r, neighbors, distances);
    LADY_BODY_COORD y(domain->getNumDim());
    for (int i = 0; i < neighbors[0].size(); ++i) {
        int cellIdx = cellCoordsMap[neighbors[0][i]];
        TracerMeshCell *cell = &tracerMeshCells(timeIdx, cellIdx);
        // calculate the tracer shape function for the cell
        (*tracer)->getBodyCoord(*domain, timeIdx, cell->getCoord(), y);
        double f = (*tracer)->getShapeFunction(timeIdx, y);
        if (f > 0.0) {
            cell->connect(*tracer, f);
            (*tracer)->connect(cell, f);
        }
    }
    if ((*tracer)->getConnectedCells().size() == 0) {
        // tracer has not connect with any cells, so connect with its host cell
        double weight = ShapeFunction::getMaxValue()/(*tracer)->getDetH(timeIdx);
        (*tracer)->getHostCell()->connect(*tracer, weight);
        (*tracer)->connect((*tracer)->getHostCell(), weight);
        if ((*tracer)->getDetH(timeIdx) < (*tracer)->getHostCell()->getVolume()) {
            // tracer is not resolved by the mesh
            (*tracer)->setBadType(Tracer::NOT_RESOLVED);
            if (numUnresolvedTracer == unresolvedTracers.size()) {
                unresolvedTracers.push_back(tracer);
            } else {
                unresolvedTracers[numUnresolvedTracer] = tracer;
            }
            numUnresolvedTracer++;
        } else {
#ifndef NDEBUG
            (*tracer)->outputNeighbors(timeIdx, *domain);
            CHECK_POINT;
#endif
        }
    } else if ((*tracer)->getBadType() == Tracer::NOT_RESOLVED) {
#ifndef NDEBUG
        (*tracer)->outputNeighbors(timeIdx, *domain);
        CHECK_POINT;
#endif
        (*tracer)->setBadType(Tracer::GOOD_SHAPE);
    }
}

void AdvectionManager::connectTracersAndMesh(const TimeLevelIndex<2> &timeIdx) {
    numUnresolvedTracer = 0;
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
        tracerMeshCells(timeIdx, i).resetConnectedTracers();
    }
    LADY_LIST<Tracer*>::iterator tracer = tracerManager.tracers.begin();
    for (; tracer != tracerManager.tracers.end(); ++tracer) {
        (*tracer)->resetConnectedCells();
        connectTracerAndMesh(timeIdx, tracer);
    }
}

void AdvectionManager::splitTracers(const TimeLevelIndex<2> &timeIdx) {
    for (int t = 0; t < numLongTracer; ++t) {
        list<Tracer*>::iterator tracer = longTracers[t];
#ifndef NDEBUG
        (*tracer)->outputNeighbors(timeIdx, *domain);
#endif
        // ---------------------------------------------------------------------
        // add two new tracers to replace the needle tracer
        Tracer *tracer1 = new Tracer(domain->getNumDim());
        Tracer *tracer2 = new Tracer(domain->getNumDim());
        const mat &H = (*tracer)->getH(timeIdx);
        const mat &U = (*tracer)->getU();
        const mat &V = (*tracer)->getV();
        vec S = (*tracer)->getS(); S[0] *= 0.5;
        // place the two tracers along major axis of the needle tracer
        static const double scale = 1/sqrt(2.0);
        LADY_BODY_COORD y(domain->getNumDim());
        y() = H*V.col(0)*scale;
        (*tracer)->getSpaceCoord(*domain, timeIdx, y, tracer1->getX(timeIdx));
        tracer1->getX(timeIdx).transformToCart(*domain);
        tracer1->getMeshIndex(timeIdx).locate(*mesh, tracer1->getX(timeIdx));
        tracer1->getH(timeIdx) = U*diagmat(S)*V.t();
        tracer1->getU() = U;
        tracer1->getS() = S;
        tracer1->getV() = V;
        tracer1->getInvH(timeIdx) = inv(tracer1->getH(timeIdx));
        tracer1->getDetH(timeIdx) = arma::prod(S);
        tracer1->updateShapeSize(*domain, timeIdx);
        tracer1->resetSkeleton(*domain, *mesh, timeIdx);
        tracerMeshCells(timeIdx, mesh->wrapIndex(tracer1->getMeshIndex(timeIdx),
                                                 CENTER)).contain(tracer1);
        tracer1->setID(tracerManager.tracers.size());
        tracerManager.tracers.push_back(tracer1);
        connectTracerAndMesh(timeIdx, --tracerManager.tracers.end());
        y() = -y();
        (*tracer)->getSpaceCoord(*domain, timeIdx, y, tracer2->getX(timeIdx));
        tracer2->getX(timeIdx).transformToCart(*domain);
        tracer2->getMeshIndex(timeIdx).locate(*mesh, tracer2->getX(timeIdx));
        tracer2->getH(timeIdx) = tracer1->getH(timeIdx);
        tracer2->getU() = U;
        tracer2->getS() = S;
        tracer2->getV() = V;
        tracer2->getInvH(timeIdx) = tracer1->getInvH(timeIdx);
        tracer2->getDetH(timeIdx) = tracer1->getDetH(timeIdx);
        tracer2->updateShapeSize(*domain, timeIdx);
        tracer2->resetSkeleton(*domain, *mesh, timeIdx);
        tracerMeshCells(timeIdx, mesh->wrapIndex(tracer2->getMeshIndex(timeIdx),
                                                 CENTER)).contain(tracer2);
        tracer2->setID(tracerManager.tracers.size());
        tracerManager.tracers.push_back(tracer2);
        connectTracerAndMesh(timeIdx, --tracerManager.tracers.end());
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            tracer1->addSpecies();
            tracer1->getSpeciesDensity(s) = (*tracer)->getSpeciesDensity(s);
            tracer1->calcSpeciesMass(timeIdx, s);
            tracer2->addSpecies();
            tracer2->getSpeciesDensity(s) = (*tracer)->getSpeciesDensity(s);
            tracer2->calcSpeciesMass(timeIdx, s);
        }
#ifndef NDEBUG
        tracer1->outputNeighbors(timeIdx, *domain);
        tracer2->outputNeighbors(timeIdx, *domain);
#endif
        REPORT_NOTICE("Long tracer " << (*tracer)->getID() <<
                      " is split into " << tracer1->getID() << " and " <<
                      tracer2->getID() << ".");
        // ---------------------------------------------------------------------
        // disconnect the needle tracer from the cells
        for (int i = 0; i < (*tracer)->getNumConnectedCell(); ++i) {
            TracerMeshCell *cell = (*tracer)->getConnectedCells()[i];
            cell->disconnect(*tracer);
        }
        (*tracer)->getHostCell()->discontain(*tracer);
        tracerManager.tracers.erase(tracer);
    }
}

void AdvectionManager::mergeTracers(const TimeLevelIndex<2> &timeIdx) {
    // merge unresolved tracers if possible
    for (int t = 0; t < numUnresolvedTracer; ++t) {
        list<Tracer*>::iterator tracer = unresolvedTracers[t];
        TracerMeshCell *cell = (*tracer)->getHostCell();
        // ---------------------------------------------------------------------
        // check density difference
        vector<Tracer*> &tracers = cell->getContainedTracers();
        vec densityDiff(tracers.size(), arma::fill::zeros);
        LADY_BODY_COORD y(domain->getNumDim());
        for (int i = 0; i < tracers.size(); ++i) {
            if (tracers[i]->getID() != (*tracer)->getID()) {
                for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                    double rho1 = (*tracer)->getSpeciesDensity(s);
                    double rho2 = tracers[i]->getSpeciesDensity(s);
                    densityDiff[i] += fabs(rho1-rho2)/rho1;
                }
            } else {
                densityDiff[i] = 1.0e33;
            }
        }
        // ---------------------------------------------------------------------
        // merge the unsolved tracer with the tracer whose density
        // difference is smallest
        arma::uword i; densityDiff.min(i);
#ifndef NDEBUG
        tracers[i]->outputNeighbors(timeIdx, *domain);
        (*tracer)->outputNeighbors(timeIdx, *domain);
#endif
        double volume1 = (*tracer)->getDetH(timeIdx);
        double volume2 = tracers[i]->getDetH(timeIdx);
        double volume = volume1+volume2;
        const mat &U = tracers[i]->getU();
        const mat &V = tracers[i]->getV();
        vec &S = tracers[i]->getS();
        S *= volume/volume2;
        tracers[i]->getH(timeIdx) = U*diagmat(S)*V.t();
        tracers[i]->getDetH(timeIdx) = volume;
        tracers[i]->getInvH(timeIdx) = inv(tracers[i]->getH(timeIdx));
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            double mass1 = (*tracer)->getSpeciesMass(s);
            double mass2 = tracers[i]->getSpeciesMass(s);
            tracers[i]->getSpeciesDensity(s) = (mass1+mass2)/volume;
            tracers[i]->calcSpeciesMass(timeIdx, s);
        }
        REPORT_NOTICE("Small tracer " << (*tracer)->getID() <<
                      " is merged into " << tracers[i]->getID() << ".");
        // ---------------------------------------------------------------------
        // disconnect the needle tracer from the cells
        for (int i = 0; i < (*tracer)->getNumConnectedCell(); ++i) {
            TracerMeshCell *cell = (*tracer)->getConnectedCells()[i];
            cell->disconnect(*tracer);
        }
        (*tracer)->getHostCell()->discontain(*tracer);
        tracerManager.tracers.erase(tracer);
    }
}
    
void AdvectionManager::handleVoidCells(const TimeLevelIndex<2> &timeIdx) {
    // NOTE: The occurance of void cells should be rare.
    for (int c = 0; c < numVoidCell; ++c) {
        TracerMeshCell *cell = voidCells[c];
        Searcher a(cellTree, NULL, cellCoords,
                   cell->getCoord().getCartCoord(), true);
        double searchRadius = 1*RAD*domain->getRadius();
        while (true) {
            mlpack::math::Range r(0.0, searchRadius);
            vector<vector<size_t> > neighbors;
            vector<vector<double> > distances;
            a.Search(r, neighbors, distances);
            if (neighbors[0].size() != 0) {
                vector<TracerMeshCell*> ngbCells;
                for (int i = 0; i < neighbors[0].size(); ++i) {
                    int cellIdx = cellCoordsMap[neighbors[0][i]];
                    TracerMeshCell *ngbCell = &tracerMeshCells(timeIdx, cellIdx);
                    // check if the cell is not a void one
                    if (find(voidCells.begin(), voidCells.end(), ngbCell)
                        == voidCells.end()) {
                        ngbCells.push_back(ngbCell);
                    }
                }
                if (ngbCells.size() == 0) {
                    searchRadius *= 2;
                    continue;
                }
                vec weights(ngbCells.size());
                for (int i = 0; i < ngbCells.size(); ++i) {
                    double d = domain->calcDistance(cell->getCoord(),
                                                    ngbCells[i]->getCoord());
                    weights[i] = 1/d;
                }
                weights /= sum(weights);
                for (int i = 0; i < ngbCells.size(); ++i) {
                    for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                        cell->getSpeciesDensity(s) += ngbCells[i]->getSpeciesDensity(s)*weights[i];
//                        cout << setw(5) << ngbCells[i]->getID();
//                        cout << setw(20) << setprecision(10) << ngbCells[i]->getSpeciesDensity(s);
//                        cout << setw(20) << setprecision(10) << weights[i];
//                        cout << setw(20) << setprecision(10) << cell->getSpeciesDensity(s) << endl;
                    }
                }
                for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                    cell->getSpeciesMass(s) = cell->getSpeciesDensity(s)*cell->getVolume();
                }
                break;
            } else {
                searchRadius *= 2;
            }
        }
    }
    
}

#define REMAP_DENSITY

void AdvectionManager::remapMeshToTracers(const TimeLevelIndex<2> &timeIdx) {
    LADY_LIST<Tracer*>::iterator tracer = tracerManager.tracers.begin();
    for (; tracer != tracerManager.tracers.end(); ++tracer) {
        (*tracer)->resetSpecies();
        const vector<TracerMeshCell*> &cells = (*tracer)->getConnectedCells();
        double totalWeight = 0;
#if defined REMAP_DENSITY
        for (int i = 0; i < (*tracer)->getNumConnectedCell(); ++i) {
            double weight = cells[i]->getRemapWeight(*tracer)*cells[i]->getVolume();
            totalWeight += weight;
            for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                (*tracer)->getSpeciesDensity(s) += cells[i]->getSpeciesDensity(s)*weight;
            }
        }
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            (*tracer)->getSpeciesDensity(s) /= totalWeight;
            (*tracer)->calcSpeciesMass(timeIdx, s);
        }
#elif defined REMAP_MASS
        for (int i = 0; i < (*tracer)->getNumConnectedCell(); ++i) {
            double weight = cells[i]->getRemapWeight(*tracer)*cells[i]->getVolume();
            totalWeight += weight;
            for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                (*tracer)->getSpeciesMass(s) += cells[i]->getSpeciesMass(s)*weight;
            }
        }
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            (*tracer)->getSpeciesMass(s) /= totalWeight;
            (*tracer)->getSpeciesDensity(s) = (*tracer)->getSpeciesMass(s)/(*tracer)->getDetH(timeIdx);
        }
#endif
//        if ((*tracer)->getID() == 0) {
//            for (int i = 0; i < (*tracer)->getNumConnectedCell(); ++i) {
//                double weight = cells[i]->getRemapWeight(*tracer)*cells[i]->getVolume()/totalWeight;
//                cout << setw(8) << cells[i]->getID();
//                cout << setw(40) << cells[i]->getSpeciesDensity(0);
//                cout << setw(40) << weight;
//                cout << setw(40) << (*tracer)->getSpeciesDensity(0) << endl;
//            }
//        }
    }
}

void AdvectionManager::remapTracersToMesh(const TimeLevelIndex<2> &timeIdx) {
    numVoidCell = 0;
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
        TracerMeshCell &cell = tracerMeshCells(timeIdx, i);
        cell.resetSpecies();
        if (cell.getNumConnectedTracer() == 0) {
//            cell.getCoord().print();
            if (numVoidCell == voidCells.size()) {
                voidCells.push_back(&cell);
            } else {
                voidCells[numVoidCell] = &cell;
            }
            numVoidCell++;
            continue;
        };
        vector<Tracer*> &tracers = cell.getConnectedTracers();
        double totalWeight = 0;
#if defined REMAP_DENSITY
        for (int j = 0; j < cell.getNumConnectedTracer(); ++j) {
            double weight = cell.getRemapWeight(tracers[j])*tracers[j]->getDetH(timeIdx);
            totalWeight += weight;
            for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                cell.getSpeciesDensity(s) += tracers[j]->getSpeciesDensity(s)*weight;
            }
        }
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            cell.getSpeciesDensity(s) /= totalWeight;
            cell.getSpeciesMass(s) = cell.getSpeciesDensity(s)*cell.getVolume();
        }
#elif defined REMAP_MASS
        for (int j = 0; j < cell.getNumConnectedTracer(); ++j) {
            double weight = cell.getRemapWeight(tracers[j])*tracers[j]->getDetH(timeIdx);
            totalWeight += weight;
            for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
                cell.getSpeciesMass(s) += tracers[j]->getSpeciesMass(s)*weight;
            }
        }
        for (int s = 0; s < tracerManager.getNumSpecies(); ++s) {
            cell.getSpeciesMass(s) /= totalWeight;
            cell.getSpeciesDensity(s) = cell.getSpeciesMass(s)/cell.getVolume();
        }
#endif
//        if (cell.getID() == 0) {
//            for (int j = 0; j < cell.getNumConnectedTracer(); ++j) {
//                double weight = cell.getRemapWeight(tracers[j])*tracers[j]->getDetH(timeIdx)/totalWeight;
//                cout << setw(8) << tracers[j]->getID();
//                cout << setw(40) << tracers[j]->getSpeciesDensity(1);
//                cout << setw(40) << weight;
//                cout << setw(40) << cell.getSpeciesDensity(1) << endl;
//            }
//        }
    }
    handleVoidCells(timeIdx);
#ifdef REMAP_DENSITY
    correctTotalMassOnMesh(timeIdx);
#endif
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
#ifndef NDEBUG
        double biasPercent = (totalMass.getLevel(timeIdx)[s]-expectedTotalMass[s])/expectedTotalMass[s];
        REPORT_NOTICE("Mass conservation bias percentage is " << std::fixed << setw(10) << setprecision(4) << biasPercent*100 << "%.");
#endif
        for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
            tracerMeshCells(timeIdx, i).getSpeciesDensity(s) *= fixer;
        }
        totalMass.getLevel(timeIdx)[s] = expectedTotalMass[s];
    }
}

}
