#include "SolidRotationTestCase.h"

namespace lady {

SolidRotationTestCase::SolidRotationTestCase() {
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
    // set parameters
    // =========================================================================
    angleSpeed = PI2/12/86400;
    U0 = domain->getRadius()*angleSpeed;
    alpha = M_PI_2;
    R = domain->getRadius()/3;
    H0 = 1000;
    // -------------------------------------------------------------------------
    REPORT_ONLINE;
}

SolidRotationTestCase::~SolidRotationTestCase() {
    delete mesh;
    delete domain;
    delete axisPole;
    delete c0;
    delete cr0;
    for (int i = 0; i < q.size(); ++i) {
        delete q[i];
    }
    REPORT_OFFLINE;
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
        double cosLat = mesh->getCosLat(CENTER, j);
        double sinLat = mesh->getSinLat(CENTER, j);
        for (int i = 0; i < mesh->getNumGrid(0, EDGE); ++i) {
            double cosLon = mesh->getCosLon(EDGE, i);
            (*V)(timeLevel, 0, i, j) = U0*(cosLat*cos(alpha)+
                                           sinLat*cosLon*sin(alpha));
        }
    }
    for (int j = 0; j < mesh->getNumGrid(1, EDGE); ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, CENTER); ++i) {
            double sinLon = mesh->getSinLon(CENTER, i);
            (*V)(timeLevel, 1, i, j) = -U0*sinLon*sin(alpha);
        }
    }
    // TEST: calculate velocity gradient tensor analytically
#define USE_ANALYTICAL_TENSOR 0
#if USE_ANALYTICAL_TENSOR == 1
    for (int j = 0; j < mesh->getNumGrid(1, CENTER); ++j) {
        double cosLat = mesh->getCosLat(CENTER, j);
        double sinLat = mesh->getSinLat(CENTER, j);
        double tanLat = mesh->getTanLat(CENTER, j);
        for (int i = 0; i < mesh->getNumGrid(0, CENTER); ++i) {
            double cosLon = mesh->getCosLon(CENTER, i);
            double sinLon = mesh->getSinLon(CENTER, i);
            double RcosLat = domain->getRadius()*cosLat;
            double u = U0*(cosLat*cos(alpha)+sinLat*cosLon*sin(alpha));
            double v = -U0*sinLon*sin(alpha);
            double dudlon = -U0*sinLat*sinLon*sin(alpha)/RcosLat;
            double dudlat = -U0*sinLat*cos(alpha)/domain->getRadius();
            double dvdlon = -U0*cosLon*sin(alpha)/RcosLat;
            double dvdlat = 0.0;
            (*T)(timeLevel, 0, 0, i, j) = dudlon;
            (*T)(timeLevel, 0, 1, i, j) = dudlat+u/domain->getRadius()*tanLat;
            (*T)(timeLevel, 1, 0, i, j) = dvdlon;
            (*T)(timeLevel, 1, 1, i, j) = dvdlat+v/domain->getRadius()*tanLat;
        }
    }
#endif
    V->applyBndCond(timeLevel, UPDATE_HALF_LEVEL);
#if USE_ANALYTICAL_TENSOR != 1
    T->calcFromVector(*V, timeLevel);
#endif
    T->applyBndCond(timeLevel, UPDATE_HALF_LEVEL);
}

}