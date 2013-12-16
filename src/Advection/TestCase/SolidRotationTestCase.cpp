#include "SolidRotationTestCase.h"

namespace lady {

SolidRotationTestCase::SolidRotationTestCase() {
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

    // set parameters
    angleSpeed = PI2/12/86400;
    U0 = domain->getRadius()*angleSpeed;
    alpha = M_PI_2;
    R = domain->getRadius()/3;
    H0 = 1000;

    REPORT_NOTICE("Solid rotation test case is online.");
}

SolidRotationTestCase::~SolidRotationTestCase() {
    delete mesh;
    delete domain;
    REPORT_NOTICE("Solid rotation test case is offline.");
}

Time SolidRotationTestCase::getStartTime() const {
    Time time;
    return time;
}

Time SolidRotationTestCase::getEndTime() const {
    Time time;
    return time+12*86400;
}

double SolidRotationTestCase::getStepSize() const {
    return 30*60;
}

void SolidRotationTestCase::advance(double time, int timeLevel) {
    for (int j = 0; j < mesh->getNumGrid(1, CENTER); ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, EDGE); ++i) {
            double lon = mesh->getGridCoord(0, EDGE, i);
            double lat = mesh->getGridCoord(1, CENTER, j);
            (*v)(timeLevel, 0, i, j) = U0*(cos(lat)*cos(alpha)+
                                           sin(lat)*cos(lon)*sin(alpha));
        }
    }
    for (int j = 0; j < mesh->getNumGrid(1, EDGE); ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, CENTER); ++i) {
            double lon = mesh->getGridCoord(0, CENTER, i);
            (*v)(timeLevel, 1, i, j) = -U0*sin(lon)*sin(alpha);
        }
    }
    v->applyBndCond(timeLevel, true);
}

}