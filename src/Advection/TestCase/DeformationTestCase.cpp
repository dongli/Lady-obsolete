#include "DeformationTestCase.h"

namespace lady {

DeformationTestCase::DeformationTestCase(SubCase subCase, InitCond initCond) {
    this->subCase = subCase;
    this->initCond = initCond;
    REPORT_ONLINE;
}

DeformationTestCase::~DeformationTestCase() {
    delete mesh;
    delete domain;
    REPORT_OFFLINE;
}

void DeformationTestCase::init(const geomtk::TimeManager &timeManager) {
    AdvectionTestCase::init(timeManager);
    period = 5.0;
    // -------------------------------------------------------------------------
    // initialize domain
    domain = new geomtk::SphereDomain(2);
    domain->setRadius(1.0);
    // -------------------------------------------------------------------------
    // initialize mesh
    mesh = new geomtk::RLLMesh(*domain);
    mesh->init(360, 181);
    // -------------------------------------------------------------------------
    // initialize velocity
    velocity.create(*mesh, true, HAS_HALF_LEVEL);
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
        for (int j = 0; j < mesh->getNumGrid(1, velocity(0).getGridType(1)); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, velocity(0).getGridType(0)); ++i) {
                double lon = mesh->getGridCoordComp(0, velocity(0).getGridType(0), i);
                double lat = mesh->getGridCoordComp(1, velocity(0).getGridType(1), j);
                velocity(0)(timeIdx, i, j) = k*pow(sin(lon*0.5), 2.0)*sin(lat*2.0)*cosT;
            }
        }
        for (int j = 0; j < mesh->getNumGrid(1, velocity(1).getGridType(1)); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, velocity(1).getGridType(0)); ++i) {
                double lon = mesh->getGridCoordComp(0, velocity(1).getGridType(0), i);
                double lat = mesh->getGridCoordComp(1, velocity(1).getGridType(1), j);
                velocity(1)(timeIdx, i, j) = k*0.5*sin(lon)*cos(lat)*cosT;
            }
        }
    } else if (subCase == CASE2) {
        k = 2.0;
        for (int j = 0; j < mesh->getNumGrid(1, velocity(0).getGridType(1)); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, velocity(0).getGridType(0)); ++i) {
                double lon = mesh->getGridCoordComp(0, velocity(0).getGridType(0), i);
                double lat = mesh->getGridCoordComp(1, velocity(0).getGridType(1), j);
                velocity(0)(timeIdx, i, j) = k*pow(sin(lon), 2.0)*sin(lat*2.0)*cosT;
            }
        }
        for (int j = 0; j < mesh->getNumGrid(1, velocity(1).getGridType(1)); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, velocity(1).getGridType(0)); ++i) {
                double lon = mesh->getGridCoordComp(0, velocity(1).getGridType(0), i);
                double lat = mesh->getGridCoordComp(1, velocity(1).getGridType(1), j);
                velocity(1)(timeIdx, i, j) = k*sin(lon*2.0)*cos(lat)*cosT;
            }
        }
    } else if (subCase == CASE3) {
        k = 1.0;
        for (int j = 0; j < mesh->getNumGrid(1, velocity(0).getGridType(1)); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, velocity(0).getGridType(0)); ++i) {
                double lon = mesh->getGridCoordComp(0, velocity(0).getGridType(0), i);
                double lat = mesh->getGridCoordComp(1, velocity(0).getGridType(1), j);
                velocity(0)(timeIdx, i, j) = -k*pow(sin(lon), 2.0)*sin(lat*2.0)*
                                        pow(cos(lat), 2.0)*cosT;
            }
        }
        for (int j = 0; j < mesh->getNumGrid(1, velocity(1).getGridType(1)); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, velocity(1).getGridType(0)); ++i) {
                double lon = mesh->getGridCoordComp(0, velocity(1).getGridType(0), i);
                double lat = mesh->getGridCoordComp(1, velocity(1).getGridType(1), j);
                velocity(1)(timeIdx, i, j) = k*0.5*sin(lon)*pow(cos(lat), 3.0)*cosT;
            }
        }
    } else if (subCase == CASE4) {
        k = 10.0*R/period;
        double c1 = PI2*time/period;
        double c2 = PI2*R/period;
        for (int j = 0; j < mesh->getNumGrid(1, velocity(0).getGridType(1)); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, velocity(0).getGridType(0)); ++i) {
                double lon = mesh->getGridCoordComp(0, velocity(0).getGridType(0), i)-c1;
                double lat = mesh->getGridCoordComp(1, velocity(0).getGridType(1), j);
                velocity(0)(timeIdx, i, j) = k*pow(sin(lon), 2.0)*sin(lat*2.0)*cosT+c2*cos(lat);
            }
        }
        for (int j = 0; j < mesh->getNumGrid(1, velocity(1).getGridType(1)); ++j) {
            for (int i = 0; i < mesh->getNumGrid(0, velocity(1).getGridType(0)); ++i) {
                double lon = mesh->getGridCoordComp(0, velocity(1).getGridType(0), i)-c1;
                double lat = mesh->getGridCoordComp(1, velocity(1).getGridType(1), j);
                velocity(1)(timeIdx, i, j) = k*sin(lon*2.0)*cos(lat)*cosT;
            }
        }
    }
    if (timeIdx.isCurrentIndex()) {
        velocity.applyBndCond(timeIdx);
    } else {
        velocity.applyBndCond(timeIdx, UPDATE_HALF_LEVEL);
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
        q0->create("", "", "", *mesh, CENTER);
        for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
            (*q0)(timeIdx, i) = 1.0;
        }
        // test tracer
        q.push_back(new LADY_SCALAR_FIELD); q1 = q.back();
        q1->create("", "", "", *mesh, CENTER);
        double b = 0.1, c = 1.0, r = 0.5;
        for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
            LADY_SPACE_COORD x(2);
            mesh->getGridCoord(i, CENTER, x);
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
