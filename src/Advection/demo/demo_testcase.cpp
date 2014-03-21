#include "lady.h"

#define USE_DEFORMATION_TEST_CASE 1
#define USE_SOLID_ROTATION_TEST_CASE 0
#define CALCULATE_SOLUTION 0

int main(int argc, const char *argv[])
{
#if USE_DEFORMATION_TEST_CASE == 1
    lady::DeformationTestCase testCase(lady::DeformationTestCase::CASE4,
                                       lady::DeformationTestCase::SLOTTED_CYLINDERS);
#elif USE_SOLID_ROTATION_TEST_CASE == 1
    lady::SolidRotationTestCase testCase;
#endif
    lady::AdvectionManager advectionManager;
    geomtk::StampString o1("tracers.deform.128x64.", ".nc");
    geomtk::TimeManager timeManager;
    geomtk::TimeLevelIndex<2> oldTimeIdx;
    // -------------------------------------------------------------------------
    // initialization
    timeManager.init(testCase.getStartTime(), testCase.getEndTime(),
                     testCase.getStepSize());

    advectionManager.init(testCase.getDomain(), testCase.getMesh(), 128, 64);

    testCase.calcInitCond(advectionManager);
    testCase.advance(timeManager.getSeconds(), oldTimeIdx);

    advectionManager.output(o1.run("%3.3d", timeManager.getNumStep()), oldTimeIdx);
    // -------------------------------------------------------------------------
    // integration loop
    while (!timeManager.isFinished()) {
        geomtk::TimeLevelIndex<2> newTimeIdx = oldTimeIdx+1;
        double time = timeManager.getSeconds()+timeManager.getStepSize();
        testCase.advance(time, newTimeIdx);
        advectionManager.advance(timeManager.getStepSize(), newTimeIdx,
                                 testCase.getVelocityField());
#if CALCULATE_SOLUTION == 1
        testCase.calcSolution(time, newTimeIdx, advectionManager);
#endif
        timeManager.advance();
        oldTimeIdx.shift();
        advectionManager.output(o1.run("%3.3d", timeManager.getNumStep()), oldTimeIdx);
    }

    return 0;
}
