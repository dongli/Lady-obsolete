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
    // Note: Use larger Pole radius to decrease deformation when crossing Poles.
    mesh->setPoleRadius(30.0*RAD);
    int numLon = 360;
    vec fullLon(numLon), halfLon(numLon);
    double dlon = 2.0*M_PI/numLon;
    for (int i = 0; i < numLon; ++i) {
        fullLon[i] = i*dlon;
        halfLon[i] = i*dlon+dlon*0.5;
    }
    mesh->setGridCoords(0, numLon, fullLon, halfLon);
    int numLat = 181;
    vec fullLat(numLat), halfLat(numLat-1);
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
    // initialize velocity
    velocity.create(*mesh, true, HAS_HALF_LEVEL);
    // -------------------------------------------------------------------------
    // set parameters
    // =========================================================================
    angleSpeed = PI2/12/86400;
    U0 = domain->getRadius()*angleSpeed;
    alpha = M_PI_2;
    // =========================================================================
    axisPole = new LADY_SPACE_COORD(2);
    c0 = new LADY_SPACE_COORD(2);
    cr0 = new LADY_SPACE_COORD(2);
    axisPole->setCoord(M_PI,  M_PI_2-alpha);
    c0->setCoord(M_PI_2, 0.0);
    domain->rotate(*axisPole, *c0, *cr0);
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

void SolidRotationTestCase::advance(double time,
                                    const TimeLevelIndex<2> &timeIdx) {
    double sinAlpha = sin(alpha), cosAlpha = cos(alpha);
    for (int j = 0; j < mesh->getNumGrid(1, velocity(0).getGridType(1)); ++j) {
        double cosLat = mesh->getCosLat(velocity(0).getGridType(1), j);
        double sinLat = mesh->getSinLat(velocity(0).getGridType(1), j);
        for (int i = 0; i < mesh->getNumGrid(0, velocity(0).getGridType(0)); ++i) {
            double cosLon = mesh->getCosLon(velocity(0).getGridType(0), i);
            velocity(0)(timeIdx, i, j) = U0*(cosLat*cosAlpha+sinLat*cosLon*sinAlpha);
        }
    }
    for (int j = 0; j < mesh->getNumGrid(1, velocity(1).getGridType(1)); ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, velocity(1).getGridType(0)); ++i) {
            double sinLon = mesh->getSinLon(velocity(1).getGridType(0), i);
            velocity(1)(timeIdx, i, j) = -U0*sinLon*sinAlpha;
        }
    }
    if (timeIdx.isCurrentIndex()) {
        velocity.applyBndCond(timeIdx);
    } else {
        velocity.applyBndCond(timeIdx, UPDATE_HALF_LEVEL);
    }
}

void SolidRotationTestCase::calcInitCond(AdvectionManager &advectionManager) {
    q.push_back(new LADY_SCALAR_FIELD);
    q.front()->create("", "", "", *mesh, CENTER);
    TimeLevelIndex<2> initTimeIdx;
    calcSolution(0, initTimeIdx, *q.front());
    AdvectionTestCase::calcInitCond(advectionManager);
}

void SolidRotationTestCase::calcSolution(double time,
                                         const TimeLevelIndex<2> &timeIdx,
                                         AdvectionManager &advectionManager) {
    calcSolution(time, timeIdx, *q.front());
    advectionManager.input(timeIdx, q);
    REPORT_NOTICE("Overwrite tracers with the true solution.");
}

void SolidRotationTestCase::calcSolution(double time,
                                         const TimeLevelIndex<2> &timeIdx,
                                         LADY_SCALAR_FIELD &q) {
    cr0->setCoordComp(0, angleSpeed*time);
    domain->rotateBack(*axisPole, *c0, *cr0);
    for (int j = 0; j < mesh->getNumGrid(1, FULL); ++j) {
        for (int i = 0; i < mesh->getNumGrid(0, FULL); ++i) {
            double lon = mesh->getGridCoordComp(0, FULL, i);
            double sinLat = mesh->getSinLat(FULL, j);
            double cosLat = mesh->getCosLat(FULL, j);
            double d = domain->calcDistance(*c0, lon, sinLat, cosLat);
            if (d < R) {
                q(timeIdx, i, j) = H0*(1+cos(M_PI*d/R))/2;
            } else {
                q(timeIdx, i, j) = 0;
            }
        }
    }
}

}
