#include "lady.h"

#define USE_DEFORMATION_TEST_CASE 0
#define USE_SOLID_ROTATION_TEST_CASE 1

int main(int argc, const char *argv[])
{
#if USE_DEFORMATION_TEST_CASE == 1
    lady::DeformationTestCase testCase(lady::DeformationTestCase::CASE4,
                                       lady::DeformationTestCase::SLOTTED_CYLINDERS);
#elif USE_SOLID_ROTATION_TEST_CASE == 1
    lady::SolidRotationTestCase testCase;
#endif
    lady::AdvectionManager advectionManager;
    geomtk::StampString o1("tracers.", ".nc"), o2("velocity.", ".nc");
    geomtk::TimeManager timeManager;
    geomtk::TimeLevelIndex<2> oldTimeIdx;
    // -------------------------------------------------------------------------
    // initialization
    timeManager.init(testCase.getStartTime(), testCase.getEndTime(),
                     testCase.getStepSize());

    advectionManager.init(testCase.getDomain(), testCase.getMesh(), 2048);

    testCase.calcInitCond(advectionManager);
    testCase.advance(timeManager.getSeconds(), oldTimeIdx);

    advectionManager.output(o1.run("%3.3d", timeManager.getNumStep()), oldTimeIdx);
    // -------------------------------------------------------------------------
    // integration loop
    while (!timeManager.isFinished()) {
        geomtk::TimeLevelIndex<2> newTimeIdx = oldTimeIdx+1;
        testCase.advance(timeManager.getSeconds()+timeManager.getStepSize(),
                         newTimeIdx);
        advectionManager.advance(timeManager.getStepSize(), newTimeIdx,
                                 testCase.getVelocityField());
        timeManager.advance();
        oldTimeIdx.shift();
        advectionManager.output(o1.run("%3.3d", timeManager.getNumStep()), oldTimeIdx);
    }

    return 0;
}
