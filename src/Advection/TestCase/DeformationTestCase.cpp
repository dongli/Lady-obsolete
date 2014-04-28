#include "DeformationTestCase.h"

namespace lady {

DeformationTestCase::DeformationTestCase(SubCase subCase, InitCond initCond) {
    this->subCase = subCase;
    this->initCond = initCond;
    period = 5.0;
    REPORT_ONLINE;
}

DeformationTestCase::~DeformationTestCase() {
    delete mesh;
    delete domain;
    REPORT_OFFLINE;
}

void DeformationTestCase::init(const geomtk::TimeManager &timeManager) {
    AdvectionTestCase::init(timeManager);
    // -------------------------------------------------------------------------
    // initialize domain
    domain = new geomtk::SphereDomain(2);
    domain->setRadius(1.0);
    // -------------------------------------------------------------------------
    // initialize mesh
    mesh = new geomtk::RLLMesh(*domain);
    mesh->init(180, 91);
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
        k = 5.0*R/period;
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
    c0.setCoord(M_PI*5.0/6.0, 0.0); c0.transformToCart(*domain);
    c1.setCoord(M_PI*7.0/6.0, 0.0); c1.transformToCart(*domain);
    // reference tracer
    q.push_back(new LADY_SCALAR_FIELD); q0 = q.back();
    q0->create("", "", "", *mesh, CENTER);
    // test tracer
    q.push_back(new LADY_SCALAR_FIELD); q1 = q.back();
    q1->create("", "", "", *mesh, CENTER);
    LADY_SPACE_COORD x(2);
    for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
        (*q0)(timeIdx, i) = 1.0;
    }
    if (initCond == COSINE_HILLS) {
        double hmax = 1, r = domain->getRadius()*0.5, b = 0.1, c = 0.9;
        for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
            mesh->getGridCoord(i, CENTER, x);
            double r0 = domain->calcDistance(x, c0);
            double r1 = domain->calcDistance(x, c1);
            if (r0 < r) {
                (*q1)(timeIdx, i) = b+c*hmax*0.5*(1+cos(M_PI*r0/r));
            } else if (r1 < r) {
                (*q1)(timeIdx, i) = b+c*hmax*0.5*(1+cos(M_PI*r1/r));
            } else {
                (*q1)(timeIdx, i) = b;
            }
        }
    } else if (initCond == SLOTTED_CYLINDERS) {
        double b = 0.1, c = 1.0, r = 0.5;
        for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
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
    } else if (initCond == GAUSSIAN_HILLS) {
        double hmax = 0.95, b = 5.0;
        for (int i = 0; i < mesh->getTotalNumGrid(CENTER); ++i) {
            mesh->getGridCoord(i, CENTER, x);
            x.transformToCart(*domain);
            vec d0 = x.getCartCoord()-c0.getCartCoord();
            vec d1 = x.getCartCoord()-c1.getCartCoord();
            (*q1)(timeIdx, i) = hmax*(exp(-b*dot(d0, d0))+exp(-b*dot(d1, d1)));
        }
    }
    AdvectionTestCase::calcInitCond(advectionManager);
}

}
