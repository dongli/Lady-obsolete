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
    return velocity;
}

void AdvectionTestCase::outputVelocity(const string &fileName,
                                       const TimeLevelIndex<2> &oldTimeIdx) const {
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