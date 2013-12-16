#include "DeformationTestCase.h"

namespace lady {

DeformationTestCase::DeformationTestCase(SubCase subCase, InitCond initCond) {
    this->subCase = subCase;
    this->initCond = initCond;
    T = 5.0;

    // initialize domain
    domain = new geomtk::SphereDomain(2);
    domain->setRadius(1.0);
    
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
    int numLat = 180;
    double fullLat[numLat], halfLat[numLat-1];
    double dlat = M_PI/(numLat-1);
    for (int j = 0; j < numLat; ++j) {
        fullLat[j] = j*dlat-M_PI_2;
    }
    for (int j = 0; j < numLat-1; ++j) {
        halfLat[j] = dlat*0.5+j*dlat-M_PI_2;
    }
    mesh->setGridCoords(1, numLat, fullLat, halfLat);

    // initialize velocity
    v = new LADY_VELOCITY_FIELD("v", "m s-1", "advection velocity",
                                *mesh, HAS_HALF_LEVEL);
    v->create(C_GRID);

    REPORT_NOTICE("Deformation test case is online.");
}

DeformationTestCase::~DeformationTestCase() {
    delete mesh;
    delete domain;

    REPORT_NOTICE("Deformation test case is offline.");
}

Time DeformationTestCase::getStartTime() const {
    Time time;
    return time;
}

Time DeformationTestCase::getEndTime() const {
    Time time = getStartTime()+1.0;
    return time;
}

double DeformationTestCase::getStepSize() const {
    return 1.0/500.0;
}

void DeformationTestCase::advance(double time, int timeLevel) {

    double cosT = cos(M_PI*time/T);
    // advance velocity
    double k;
    if (subCase == CASE1) {
        k = 2.4;
        for (int j = 0; j < mesh->getNumGrid(1, geomtk::CENTER); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, geomtk::EDGE); ++i) {
                double lon = mesh->getGridCoord(0, geomtk::EDGE, i);
                double lat = mesh->getGridCoord(1, geomtk::CENTER, j);
                (*v)(0, 0, i, j) = k*pow(sin(lon*0.5), 2.0)*sin(lat*2.0)*cosT;
            }
        }
        for (int j = 0; j < mesh->getNumGrid(1, geomtk::EDGE); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, geomtk::CENTER); ++i) {
                double lon = mesh->getGridCoord(0, geomtk::CENTER, i);
                double lat = mesh->getGridCoord(1, geomtk::EDGE, j);
                (*v)(0, 1, i, j) = k*0.5*sin(lon)*cos(lat)*cosT;
            }
        }
    } else if (subCase == CASE2) {
        k = 2.0;
        for (int j = 0; j < mesh->getNumGrid(1, geomtk::CENTER); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, geomtk::EDGE); ++i) {
                double lon = mesh->getGridCoord(0, geomtk::EDGE, i);
                double lat = mesh->getGridCoord(1, geomtk::CENTER, j);
                (*v)(0, 0, i, j) = k*pow(sin(lon), 2.0)*sin(lat*2.0)*cosT;
            }
        }
        for (int j = 0; j < mesh->getNumGrid(1, geomtk::EDGE); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, geomtk::CENTER); ++i) {
                double lon = mesh->getGridCoord(0, geomtk::CENTER, i);
                double lat = mesh->getGridCoord(1, geomtk::EDGE, j);
                (*v)(0, 1, i, j) = k*sin(lon*2.0)*cos(lat)*cosT;
            }
        }
    } else if (subCase == CASE3) {
        k = 1.0;
        for (int j = 0; j < mesh->getNumGrid(1, geomtk::CENTER); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, geomtk::EDGE); ++i) {
                double lon = mesh->getGridCoord(0, geomtk::EDGE, i);
                double lat = mesh->getGridCoord(1, geomtk::CENTER, j);
                (*v)(0, 0, i, j) = -k*pow(sin(lon), 2.0)*sin(lat*2.0)*pow(cos(lat), 2.0)*cosT;
            }
        }
        for (int j = 0; j < mesh->getNumGrid(1, geomtk::EDGE); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, geomtk::CENTER); ++i) {
                double lon = mesh->getGridCoord(0, geomtk::CENTER, i);
                double lat = mesh->getGridCoord(1, geomtk::EDGE, j);
                (*v)(0, 1, i, j) = k*0.5*sin(lon)*pow(cos(lat), 3.0)*cosT;
            }
        }
    } else if (subCase == CASE4) {
        k = 10.0*domain->getRadius()/T;
        double c1 = PI2*time/T;
        double c2 = PI2*domain->getRadius()/T;
        for (int j = 0; j < mesh->getNumGrid(1, geomtk::CENTER); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, geomtk::EDGE); ++i) {
                double lon = mesh->getGridCoord(0, geomtk::EDGE, i);
                double lat = mesh->getGridCoord(1, geomtk::CENTER, j);
                (*v)(0, 0, i, j) = k*pow(sin(lon), 2.0)*sin(lat*2.0)*cosT+c2*cos(lat);
            }
        }
        for (int j = 0; j < mesh->getNumGrid(1, geomtk::EDGE); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, geomtk::CENTER); ++i) {
                double lon = mesh->getGridCoord(0, geomtk::CENTER, i)-c1;
                double lat = mesh->getGridCoord(1, geomtk::EDGE, j);
                (*v)(0, 1, i, j) = k*sin(lon*2.0)*cos(lat)*cosT;
            }
        }
    }
    v->applyBndCond(timeLevel, UPDATE_HALF_LEVEL);
}

}