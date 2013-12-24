#include "lady.h"

using geomtk::TimeManager;
using geomtk::StampString;
using lady::DeformationTestCase;
using lady::SolidRotationTestCase;
using lady::AdvectionManager;

#define USE_DEFORMATION_TEST_CASE 1
#define USE_SOLID_ROTATION_TEST_CASE 0

int main(int argc, const char *argv[])
{
    TimeManager timeManager;
#if USE_DEFORMATION_TEST_CASE == 1
    DeformationTestCase testCase(DeformationTestCase::CASE4,
                                 DeformationTestCase::SLOTTED_CYLINDERS);
#elif USE_SOLID_ROTATION_TEST_CASE == 1
    SolidRotationTestCase testCase;
#endif
    AdvectionManager advectionManager;
    StampString o1("tracers.", ".nc"), o2("velocity.", ".nc");

    timeManager.init(testCase.getStartTime(), testCase.getEndTime(),
                     testCase.getStepSize());

    advectionManager.init(testCase.getDomain(), testCase.getMesh(), 128);

    testCase.advance(timeManager.getSeconds(), timeManager.getOldLevel());

    advectionManager.output(o1.run("%3.3d", timeManager.getNumStep()),
                            timeManager.getOldLevel());
    testCase.outputVelocity(o2.run("%3.3d", timeManager.getNumStep()),
                            timeManager.getOldLevel());

    while (!timeManager.isFinished()) {
        testCase.advance(timeManager.getSeconds()+timeManager.getStepSize(),
                         timeManager.getNewLevel());
        advectionManager.advance(timeManager.getStepSize(),
                                 timeManager.getOldLevel(),
                                 timeManager.getHalfLevel(),
                                 timeManager.getNewLevel(),
                                 testCase.getVelocityField(),
                                 testCase.getTensorField());
        timeManager.advance();
        advectionManager.output(o1.run("%3.3d", timeManager.getNumStep()),
                                timeManager.getOldLevel());
        testCase.outputVelocity(o2.run("%3.3d", timeManager.getNumStep()),
                                timeManager.getOldLevel());
    }

    return 0;
}
