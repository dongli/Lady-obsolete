#include "AdvectionTestCase.h"

namespace lady {

AdvectionTestCase::AdvectionTestCase() {
}

AdvectionTestCase::~AdvectionTestCase() {
    for (int i = 0; i < q.size(); ++i) {
        delete q[i];
    }
}

const LADY_DOMAIN& AdvectionTestCase::getDomain() const {
    return *domain;
}

const LADY_MESH& AdvectionTestCase::getMesh() const {
    return *mesh;
}

const LADY_VELOCITY_FIELD& AdvectionTestCase::getVelocityField() const {
    return *V;
}

void AdvectionTestCase::outputVelocity(const string &fileName,
                                       const TimeLevelIndex<2> &oldTimeIdx) const {
    if (dynamic_cast<geomtk::RLLMesh*>(mesh)) {
        int ncId, lonDimId, latDimId, levDimId, dimDimId;
        int lonVarId, latVarId, levVarId;
        int dimIds[5];
        int uVarId, vVarId, wVarId;
        vec lon, lat, lev;
        char s[100];

        if (nc_create(fileName.c_str(), NC_CLOBBER, &ncId)
            != NC_NOERR) {
            REPORT_ERROR("Failed to create " << fileName << "!");
        }
        // ---------------------------------------------------------------------
        // define dimensions
        // =====================================================================
        // longitude dimension
        if (nc_def_dim(ncId, "lon", mesh->getNumGrid(0, CENTER), &lonDimId)
            != NC_NOERR) {
            REPORT_ERROR("Failed to define dimension lon!");
        }
        if (nc_def_var(ncId, "lon", NC_DOUBLE, 1, &lonDimId, &lonVarId)
            != NC_NOERR) {
            REPORT_ERROR("Failed to define coordinate variable lon!");
        }
        sprintf(s, "longitude");
        if (nc_put_att(ncId, lonVarId, "long_name", NC_CHAR, strlen(s), s)
            != NC_NOERR) {
            REPORT_ERROR("Failed to put attribute to variable lon!");
        }
        sprintf(s, "degree_east");
        if (nc_put_att(ncId, lonVarId, "units", NC_CHAR, strlen(s), s)
            != NC_NOERR) {
            REPORT_ERROR("Failed to put attribute to variable lon!");
        }
        // =====================================================================
        // latitude dimension
        if (nc_def_dim(ncId, "lat", mesh->getNumGrid(1, CENTER), &latDimId)
            != NC_NOERR) {
            REPORT_ERROR("Failed to define dimension lat!");
        }
        if (nc_def_var(ncId, "lat", NC_DOUBLE, 1, &latDimId, &latVarId)
            != NC_NOERR) {
            REPORT_ERROR("Failed to define coordinate variable lat!");
        }
        sprintf(s, "latitude");
        if (nc_put_att(ncId, latVarId, "long_name", NC_CHAR, strlen(s), s)
            != NC_NOERR) {
            REPORT_ERROR("Failed to put attribute to variable lat!");
        }
        sprintf(s, "degree_north");
        if (nc_put_att(ncId, latVarId, "units", NC_CHAR, strlen(s), s)
            != NC_NOERR) {
            REPORT_ERROR("Failed to put attribute to variable lat!");
        }
        // =====================================================================
        // level dimension
        if (domain->getNumDim() == 3) {
            if (nc_def_dim(ncId, "lev", mesh->getNumGrid(2, CENTER), &levDimId)
                != NC_NOERR) {
                REPORT_ERROR("Failed to define dimension lev!");
            }
            if (nc_def_var(ncId, "lev", NC_DOUBLE, 1, &levDimId, &levVarId)
                != NC_NOERR) {
                REPORT_ERROR("Failed to define coordinate variable lev!");
            }
            sprintf(s, "level");
            if (nc_put_att(ncId, levVarId, "long_name", NC_CHAR, strlen(s), s)
                != NC_NOERR) {
                REPORT_ERROR("Failed to put attribute to variable lev");
            }
        }
        // =====================================================================
        // dimension dimension
        if (nc_def_dim(ncId, "dim", domain->getNumDim(), &dimDimId)
            != NC_NOERR) {
            REPORT_ERROR("Failed to define dimension dim!");
        }
        // ---------------------------------------------------------------------
        // define variables
        if (domain->getNumDim() == 2) {
            dimIds[1] = lonDimId;
            dimIds[0] = latDimId;
        } else if (domain->getNumDim() == 3) {
            dimIds[2] = lonDimId;
            dimIds[1] = latDimId;
            dimIds[0] = levDimId;
        }
        // =====================================================================
        if (nc_def_var(ncId, "u", NC_DOUBLE, domain->getNumDim(), dimIds,
            &uVarId) != NC_NOERR) {
            REPORT_ERROR("Failed to define variable u!");
        }
        // =====================================================================
        if (nc_def_var(ncId, "v", NC_DOUBLE, domain->getNumDim(), dimIds,
            &vVarId) != NC_NOERR) {
            REPORT_ERROR("Failed to define variable v!");
        }
        // =====================================================================
        if (domain->getNumDim() == 3) {
            if (nc_def_var(ncId, "w", NC_DOUBLE, domain->getNumDim(), dimIds,
                &wVarId) != NC_NOERR) {
                REPORT_ERROR("Failed to define variable v!");
            }
        }
        // =====================================================================
        if (domain->getNumDim() == 2) {
            dimIds[3] = dimDimId;
            dimIds[2] = dimDimId;
            dimIds[1] = lonDimId;
            dimIds[0] = latDimId;
        } else if (domain->getNumDim() == 3) {
            dimIds[4] = dimDimId;
            dimIds[3] = dimDimId;
            dimIds[2] = lonDimId;
            dimIds[1] = latDimId;
            dimIds[0] = levDimId;
        }
        // ---------------------------------------------------------------------
        if (nc_enddef(ncId)) {
            REPORT_ERROR("Failed to end define!");
        }
        // ---------------------------------------------------------------------
        // put variables
        // =====================================================================
        lon = mesh->getGridCoords(0, CENTER)/RAD;
        lat = mesh->getGridCoords(1, CENTER)/RAD;
        if (nc_put_var(ncId, lonVarId, lon.memptr()) != NC_NOERR) {
            REPORT_ERROR("Failed to put coordinate variable lon!");
        }
        if (nc_put_var(ncId, latVarId, lat.memptr()) != NC_NOERR) {
            REPORT_ERROR("Failed to put coordinate variable lat!");
        }
        // =====================================================================
        if (domain->getNumDim() == 2) {
            mat u, v;
            u.reshape(mesh->getNumGrid(0, CENTER), mesh->getNumGrid(1, CENTER), 1);
            v.reshape(mesh->getNumGrid(0, CENTER), mesh->getNumGrid(1, CENTER), 1);
            V->convert(A_GRID, oldTimeIdx, u, v);
            if (nc_put_var(ncId, uVarId, u.memptr()) != NC_NOERR) {
                REPORT_ERROR("Failed to put variable u!");
            }
            if (nc_put_var(ncId, vVarId, v.memptr()) != NC_NOERR) {
                REPORT_ERROR("Failed to put variable v!");
            }
        } else if (domain->getNumDim() == 3) {
            cube u, v, w;
            u.reshape(mesh->getNumGrid(0, CENTER), mesh->getNumGrid(1, CENTER),
                      mesh->getNumGrid(2, CENTER));
            v.reshape(mesh->getNumGrid(0, CENTER), mesh->getNumGrid(1, CENTER),
                      mesh->getNumGrid(2, CENTER));
            w.reshape(mesh->getNumGrid(0, CENTER), mesh->getNumGrid(1, CENTER),
                      mesh->getNumGrid(2, CENTER));
            V->convert(A_GRID, oldTimeIdx, u, v, w);
            if (nc_put_var(ncId, uVarId, u.memptr()) != NC_NOERR) {
                REPORT_ERROR("Failed to put variable u!");
            }
            if (nc_put_var(ncId, vVarId, v.memptr()) != NC_NOERR) {
                REPORT_ERROR("Failed to put variable v!");
            }
            if (nc_put_var(ncId, wVarId, w.memptr()) != NC_NOERR) {
                REPORT_ERROR("Failed to put variable w!");
            }
        }
        // ---------------------------------------------------------------------
        if (nc_close(ncId)) {
            REPORT_ERROR("Failed to close file " << fileName << "!");
        }
    }
}

void AdvectionTestCase::calcInitCond(AdvectionManager &advectionManager) {
    assert(q.size() > 0);
    geomtk::StampString name("q", ""), brief("test tracer ", "");
    for (int i = 0; i < q.size(); ++i) {
        advectionManager.registerTracer(name.run("%d", i), "N/A",
                                        brief.run("%d", i));
    }
    TimeLevelIndex<2> initTimeIdx;
    advectionManager.input(initTimeIdx, q);
}

void AdvectionTestCase::calcSolution(double time,
                                     const TimeLevelIndex<2> &timeIdx,
                                     AdvectionManager &advectionManager) {
}

void AdvectionTestCase::calcSolution(double time,
                                     const TimeLevelIndex<2> &timeIdx,
                                     LADY_SCALAR_FIELD &q) {
    REPORT_ERROR("calcSolution is not available!");
}

}