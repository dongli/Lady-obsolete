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
    // -------------------------------------------------------------------------
    // initialize velocity and its gradient tensor
    V = new LADY_VELOCITY_FIELD("v", "m s-1", "advection velocity",
                                *mesh, HAS_HALF_LEVEL);
    V->create(C_GRID);
    T = new LADY_TENSOR_FIELD("t", "s-1", "velocity gradient tensor",
                              *mesh, HAS_HALF_LEVEL);
    T->create();
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

void DeformationTestCase::advance(double time, int timeLevel) {
    double cosT = cos(M_PI*time/period);
    double k, R = domain->getRadius();
    // advance velocity
    if (subCase == CASE1) {
        k = 2.4;
        for (int j = 0; j < mesh->getNumGrid(1, CENTER); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, EDGE); ++i) {
                double lon = mesh->getGridCoord(0, EDGE, i);
                double lat = mesh->getGridCoord(1, CENTER, j);
                (*V)(timeLevel, 0, i, j) =
                    k*pow(sin(lon*0.5), 2.0)*sin(lat*2.0)*cosT;
            }
        }
        for (int j = 0; j < mesh->getNumGrid(1, EDGE); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, CENTER); ++i) {
                double lon = mesh->getGridCoord(0, CENTER, i);
                double lat = mesh->getGridCoord(1, EDGE, j);
                (*V)(timeLevel, 1, i, j) = k*0.5*sin(lon)*cos(lat)*cosT;
            }
        }
    } else if (subCase == CASE2) {
        k = 2.0;
        for (int j = 0; j < mesh->getNumGrid(1, CENTER); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, EDGE); ++i) {
                double lon = mesh->getGridCoord(0, EDGE, i);
                double lat = mesh->getGridCoord(1, CENTER, j);
                (*V)(timeLevel, 0, i, j) =
                    k*pow(sin(lon), 2.0)*sin(lat*2.0)*cosT;
            }
        }
        for (int j = 0; j < mesh->getNumGrid(1, EDGE); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, CENTER); ++i) {
                double lon = mesh->getGridCoord(0, CENTER, i);
                double lat = mesh->getGridCoord(1, EDGE, j);
                (*V)(timeLevel, 1, i, j) = k*sin(lon*2.0)*cos(lat)*cosT;
            }
        }
    } else if (subCase == CASE3) {
        k = 1.0;
        for (int j = 0; j < mesh->getNumGrid(1, CENTER); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, EDGE); ++i) {
                double lon = mesh->getGridCoord(0, EDGE, i);
                double lat = mesh->getGridCoord(1, CENTER, j);
                (*V)(timeLevel, 0, i, j) =
                    -k*pow(sin(lon), 2.0)*sin(lat*2.0)*pow(cos(lat), 2.0)*cosT;
            }
        }
        for (int j = 0; j < mesh->getNumGrid(1, EDGE); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, CENTER); ++i) {
                double lon = mesh->getGridCoord(0, CENTER, i);
                double lat = mesh->getGridCoord(1, EDGE, j);
                (*V)(timeLevel, 1, i, j) =
                    k*0.5*sin(lon)*pow(cos(lat), 3.0)*cosT;
            }
        }
    } else if (subCase == CASE4) {
        k = 10.0*R/period;
        double c1 = PI2*time/period;
        double c2 = PI2*R/period;
        for (int j = 0; j < mesh->getNumGrid(1, CENTER); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, EDGE); ++i) {
                double lon = mesh->getGridCoord(0, EDGE, i)-c1;
                double lat = mesh->getGridCoord(1, CENTER, j);
                (*V)(timeLevel, 0, i, j) =
                    k*pow(sin(lon), 2.0)*sin(lat*2.0)*cosT+c2*cos(lat);
            }
        }
        for (int j = 0; j < mesh->getNumGrid(1, EDGE); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, CENTER); ++i) {
                double lon = mesh->getGridCoord(0, CENTER, i)-c1;
                double lat = mesh->getGridCoord(1, EDGE, j);
                (*V)(timeLevel, 1, i, j) = k*sin(lon*2.0)*cos(lat)*cosT;
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
    V->applyBndCond(timeLevel, UPDATE_HALF_LEVEL);
#if USE_ANALYTICAL_TENSOR != 1
    T->calcFromVector(*V, timeLevel);
#endif
    T->applyBndCond(timeLevel, UPDATE_HALF_LEVEL);
}

}