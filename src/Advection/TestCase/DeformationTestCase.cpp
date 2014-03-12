#include "DeformationTestCase.h"

namespace lady {

DeformationTestCase::DeformationTestCase(SubCase subCase, InitCond initCond) {
    this->subCase = subCase;
    this->initCond = initCond;
    period = 5.0;
    // -------------------------------------------------------------------------
    // initialize domain
    domain = new geomtk::SphereDomain(2);
    domain->setRadius(1.0);
    // -------------------------------------------------------------------------
    // initialize mesh
    mesh = new geomtk::RLLMesh(*domain);
    int numLon = 360;
    double fullLon[numLon], halfLon[numLon];
    double dlon = 2.0*M_PI/numLon;
    for (int i = 0; i < numLon; ++i) {
        fullLon[i] = i*dlon;
        halfLon[i] = i*dlon+dlon*0.5;
    }
    mesh->setGridCoords(0, numLon, fullLon, halfLon);
    int numLat = 181;
    double fullLat[numLat], halfLat[numLat-1];
    // NOTE: Since the velocity interpolation within polar cap is inaccurate in
    //       deformational flow, we set the second and last second latitude very
    //       near to Poles by the parameter 'lat0' (distance from Poles).
    double lat0 = 0.1*RAD;
    double dlat = (M_PI-lat0*2)/(numLat-2-1);
    for (int j = 1; j < numLat-1; ++j) {
        fullLat[j] = (j-1)*dlat-M_PI_2+lat0;
    }
    fullLat[0] = -M_PI_2;
    fullLat[numLat-1] = M_PI_2;
    for (int j = 1; j < numLat-2; ++j) {
        halfLat[j] = dlat*0.5+(j-1)*dlat-M_PI_2+lat0;
    }
    halfLat[0] = -M_PI_2+lat0*0.5;
    halfLat[numLat-2] = M_PI_2-lat0*0.5;
    mesh->setGridCoords(1, numLat, fullLat, halfLat);
    mesh->setCellVolumes();
    // -------------------------------------------------------------------------
    // initialize velocity and its gradient tensor
    V.create("v", "m s-1", "advection velocity", *mesh, _2D, C_GRID, HAS_HALF_LEVEL);
    // -------------------------------------------------------------------------
    REPORT_ONLINE;
}

DeformationTestCase::~DeformationTestCase() {
    delete mesh;
    delete domain;
    REPORT_OFFLINE;
}

Time DeformationTestCase::getStartTime() const {
    Time time;
    return time;
}

Time DeformationTestCase::getEndTime() const {
    Time time = getStartTime()+period;
    return time;
}

double DeformationTestCase::getStepSize() const {
    return period/600.0;
}

void DeformationTestCase::advance(double time,
                                  const TimeLevelIndex<2> &timeIdx) {
    double cosT = cos(M_PI*time/period);
    double k, R = domain->getRadius();
    // advance velocity
    if (subCase == CASE1) {
        k = 2.4;
        for (int j = 0; j < mesh->getNumGrid(1, CENTER); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, EDGE); ++i) {
                double lon = mesh->getGridCoordComp(0, EDGE, i);
                double lat = mesh->getGridCoordComp(1, CENTER, j);
                V(0, timeIdx, i, j) = k*pow(sin(lon*0.5), 2.0)*sin(lat*2.0)*cosT;
            }
        }
        for (int j = 0; j < mesh->getNumGrid(1, EDGE); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, CENTER); ++i) {
                double lon = mesh->getGridCoordComp(0, CENTER, i);
                double lat = mesh->getGridCoordComp(1, EDGE, j);
                V(1, timeIdx, i, j) = k*0.5*sin(lon)*cos(lat)*cosT;
            }
        }
    } else if (subCase == CASE2) {
        k = 2.0;
        for (int j = 0; j < mesh->getNumGrid(1, CENTER); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, EDGE); ++i) {
                double lon = mesh->getGridCoordComp(0, EDGE, i);
                double lat = mesh->getGridCoordComp(1, CENTER, j);
                V(0, timeIdx, i, j) = k*pow(sin(lon), 2.0)*sin(lat*2.0)*cosT;
            }
        }
        for (int j = 0; j < mesh->getNumGrid(1, EDGE); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, CENTER); ++i) {
                double lon = mesh->getGridCoordComp(0, CENTER, i);
                double lat = mesh->getGridCoordComp(1, EDGE, j);
                V(1, timeIdx, i, j) = k*sin(lon*2.0)*cos(lat)*cosT;
            }
        }
    } else if (subCase == CASE3) {
        k = 1.0;
        for (int j = 0; j < mesh->getNumGrid(1, CENTER); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, EDGE); ++i) {
                double lon = mesh->getGridCoordComp(0, EDGE, i);
                double lat = mesh->getGridCoordComp(1, CENTER, j);
                V(0, timeIdx, i, j) = -k*pow(sin(lon), 2.0)*sin(lat*2.0)*
                                        pow(cos(lat), 2.0)*cosT;
            }
        }
        for (int j = 0; j < mesh->getNumGrid(1, EDGE); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, CENTER); ++i) {
                double lon = mesh->getGridCoordComp(0, CENTER, i);
                double lat = mesh->getGridCoordComp(1, EDGE, j);
                V(1, timeIdx, i, j) = k*0.5*sin(lon)*pow(cos(lat), 3.0)*cosT;
            }
        }
    } else if (subCase == CASE4) {
        k = 10.0*R/period;
        double c1 = PI2*time/period;
        double c2 = PI2*R/period;
        for (int j = 0; j < mesh->getNumGrid(1, CENTER); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, EDGE); ++i) {
                double lon = mesh->getGridCoordComp(0, EDGE, i)-c1;
                double lat = mesh->getGridCoordComp(1, CENTER, j);
                V(0, timeIdx, i, j) = k*pow(sin(lon), 2.0)*sin(lat*2.0)*cosT+c2*cos(lat);
            }
        }
        for (int j = 0; j < mesh->getNumGrid(1, EDGE); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, CENTER); ++i) {
                double lon = mesh->getGridCoordComp(0, CENTER, i)-c1;
                double lat = mesh->getGridCoordComp(1, EDGE, j);
                V(1, timeIdx, i, j) = k*sin(lon*2.0)*cos(lat)*cosT;
            }
        }
        // TEST: calculate velocity gradient tensor analytically
#define USE_ANALYTICAL_TENSOR 0
#if USE_ANALYTICAL_TENSOR == 1
        for (int j = 0; j < mesh->getNumGrid(1, CENTER); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, CENTER); ++i) {
                double lon = mesh->getGridCoord(0, CENTER, i)-c1;
                double lat = mesh->getGridCoord(1, CENTER, j);
                double RcosLat = R*cos(lat);
                double tanLat = tan(lat);
                double u = k*pow(sin(lon), 2.0)*sin(lat*2.0)*cosT+c2*cos(lat);
                double v = k*sin(lon*2.0)*cos(lat)*cosT;
                double dudlon = k*2.0*sin(lon)*cos(lon)*sin(lat*2.0)*cosT/RcosLat;
                double dudlat = (k*pow(sin(lon), 2.0)*cos(lat*2.0)*2.0*cosT-c2*sin(lat))/R;
                double dvdlon = k*cos(lon*2.0)*2.0*cos(lat)*cosT/RcosLat;
                double dvdlat = -k*sin(lon*2.0)*sin(lat)*cosT/R;
                (*T)(timeLevel, 0, 0, i, j) = dudlon;
                (*T)(timeLevel, 0, 1, i, j) = dudlat+u/R*tanLat;
                (*T)(timeLevel, 1, 0, i, j) = dvdlon;
                (*T)(timeLevel, 1, 1, i, j) = dvdlat+v/R*tanLat;
            }
        }
#endif
    }
    if (timeIdx.isCurrentIndex()) {
        V.applyBndCond(timeIdx);
    } else {
        V.applyBndCond(timeIdx, UPDATE_HALF_LEVEL);
    }
}

void DeformationTestCase::calcInitCond(AdvectionManager &advectionManager) {
    LADY_SCALAR_FIELD *q0, *q1;
    TimeLevelIndex<2> timeIdx;
    LADY_SPACE_COORD c0(2), c1(2);
    c0.setCoord(M_PI*5.0/6.0, 0.0);
    c1.setCoord(M_PI*7.0/6.0, 0.0);
    if (initCond == COSINE_HILL) {
        
    } else if (initCond == SLOTTED_CYLINDERS) {
        // reference tracer
        q.push_back(new LADY_SCALAR_FIELD); q0 = q.back();
        q0->create("", "", "", *mesh, ScalarField, 2, A_GRID);
        for (int i = 0; i < mesh->getTotalNumGrid(A_GRID); ++i) {
            (*q0)(timeIdx, i) = 1.0;
        }
        // test tracer
        q.push_back(new LADY_SCALAR_FIELD); q1 = q.back();
        q1->create("", "", "", *mesh, ScalarField, 2, A_GRID);
        double b = 0.1, c = 1.0, r = 0.5;
        for (int i = 0; i < mesh->getTotalNumGrid(A_GRID); ++i) {
            LADY_SPACE_COORD x(2);
            mesh->getGridCoord(i, x, A_GRID);
            double r0 = domain->calcDistance(x, c0);
            double r1 = domain->calcDistance(x, c1);
            if ((r0 <= r && fabs(x(0)-c0(0)) >= r/6.0) ||
                (r1 <= r && fabs(x(0)-c1(0)) >= r/6.0))
                (*q1)(timeIdx, i) = c;
            else if (r0 <= r && fabs(x(0)-c0(0)) < r/6.0 &&
                     x(1)-c0(1) < -5.0/12.0*r)
                (*q1)(timeIdx, i) = c;
            else if (r1 <= r && fabs(x(0)-c1(0)) < r/6.0 &&
                     x(1)-c1(1) > 5.0/12.0*r)
                (*q1)(timeIdx, i) = c;
            else
                (*q1)(timeIdx, i) = b;
        }
    } else if (initCond == GAUSSIAN_HILL) {
        REPORT_ERROR("Under construction!");
    }
    AdvectionTestCase::calcInitCond(advectionManager);
}

}
