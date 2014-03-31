#include "BarotropicTestCase.h"

namespace lady {

BarotropicTestCase::BarotropicTestCase() {
    REPORT_ONLINE;
}

BarotropicTestCase::~BarotropicTestCase() {
    REPORT_OFFLINE;
}

void BarotropicTestCase::init(const geomtk::TimeManager &timeManager) {
    AdvectionTestCase::init(timeManager);
    model.init(80, 41);
    velocity.create(model.getMesh(), false, HAS_HALF_LEVEL);

    io.init(timeManager);
    fileIdx = io.registerOutputFile(model.getMesh(), "barotropic-output",
                                    geomtk::IOFrequencyUnit::STEPS, 1);
    io.file(fileIdx).registerOutputField<double, 2, barotropic_model::FULL_DIMENSION>(1, &model.getGeopotentialDepth());
}

Time BarotropicTestCase::getStartTime() const {
    Time time;
    return time;
}

Time BarotropicTestCase::getEndTime() const {
    Time time;
    return time+1*TimeUnit::DAYS;
}

double BarotropicTestCase::getStepSize() const {
    return 2*TimeUnit::MINUTES;
}

void BarotropicTestCase::calcInitCond(AdvectionManager &advectionManager) {
    // set initial condition for barotropic model
    LADY_SPACE_COORD x(2);
    x.setCoord(180*RAD, 30*RAD);
    testCase.addPeak(x, 1000*barotropic_model::G, model.getDomain().getRadius()/3);
    testCase.calcInitCond(model);
    // set initial condition for tracer
    const LADY_SCALAR_FIELD &gd = model.getGeopotentialDepth();
    q.push_back(new LADY_SCALAR_FIELD);
    q.front()->create("", "", "", model.getMesh(), CENTER);
    TimeLevelIndex<2> initTimeIdx;
    for (int j = 0; j < model.getMesh().getNumGrid(1, FULL); ++j) {
        for (int i = 0; i < model.getMesh().getNumGrid(0, FULL); ++i) {
            (*q.front())(initTimeIdx, i, j) = gd(initTimeIdx, i, j);
        }
    }
    AdvectionTestCase::calcInitCond(advectionManager);
}

void BarotropicTestCase::advance(double time,
                                 const TimeLevelIndex<2> &timeIdx) {
    if (timeIdx.isCurrentIndex()) {
        model.integrate(timeIdx, getStepSize());
    } else {
        model.integrate(timeIdx-1, getStepSize());
    }
    io.create(fileIdx);
    io.output<double, 2>(fileIdx, timeIdx, 1, &model.getGeopotentialDepth());
    io.close(fileIdx);
    for (int j = 0; j < model.getMesh().getNumGrid(1, velocity(0).getGridType(1)); ++j) {
        for (int i = 0; i < model.getMesh().getNumGrid(0, velocity(0).getGridType(0)); ++i) {
            velocity(0)(timeIdx, i, j) = model.getZonalWind()(timeIdx, i, j);
        }
    }
    for (int j = 0; j < model.getMesh().getNumGrid(1, velocity(1).getGridType(1)); ++j) {
        for (int i = 0; i < model.getMesh().getNumGrid(0, velocity(1).getGridType(0)); ++i) {
            velocity(1)(timeIdx, i, j) = model.getMeridionalWind()(timeIdx, i, j);
        }
    }
    if (timeIdx.isCurrentIndex()) {
        velocity.applyBndCond(timeIdx);
    } else {
        velocity.applyBndCond(timeIdx, UPDATE_HALF_LEVEL);
    }
}
    
}