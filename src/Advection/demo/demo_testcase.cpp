#include "lady.h"

#define USE_DEFORMATION_TEST_CASE 0
#define USE_SOLID_ROTATION_TEST_CASE 0
#define USE_BAROTROPIC_TEST_CASE 1
#define CALCULATE_SOLUTION 0

int main(int argc, const char *argv[])
{
#if USE_DEFORMATION_TEST_CASE == 1
    lady::DeformationTestCase testCase(lady::DeformationTestCase::CASE4,
                                       lady::DeformationTestCase::SLOTTED_CYLINDERS);
    geomtk::StampString o1("tracers.deform.240x120.", ".nc");
#elif USE_SOLID_ROTATION_TEST_CASE == 1
    lady::SolidRotationTestCase testCase;
#if CALCULATE_SOLUTION == 1
    geomtk::StampString o1("tracers.rotation.true.240x120.", ".nc");
#else
    geomtk::StampString o1("tracers.rotation.480x240.", ".nc");
#endif
#elif USE_BAROTROPIC_TEST_CASE == 1
    lady::BarotropicTestCase testCase;
    geomtk::StampString o1("tracers.barotropic.80x40.", ".nc");
#endif
    lady::AdvectionManager advectionManager;
    geomtk::TimeManager timeManager;
    geomtk::TimeLevelIndex<2> oldTimeIdx;
    // -------------------------------------------------------------------------
    // initialization
    timeManager.init(testCase.getStartTime(), testCase.getEndTime(),
                     testCase.getStepSize());
    testCase.init(timeManager);
    advectionManager.init(testCase.getDomain(), testCase.getMesh(), 80, 40);

    testCase.calcInitCond(advectionManager);
    advectionManager.output(o1.run("%4.4d", timeManager.getNumStep()), oldTimeIdx);

    testCase.advance(timeManager.getSeconds(), oldTimeIdx);
    // -------------------------------------------------------------------------
    // integration loop
    while (!timeManager.isFinished()) {
        geomtk::TimeLevelIndex<2> newTimeIdx = oldTimeIdx+1;
        double time = timeManager.getSeconds()+timeManager.getStepSize();
        testCase.advance(time, newTimeIdx);
        advectionManager.advance(timeManager.getStepSize(), newTimeIdx,
                                 testCase.getVelocityField());
#if CALCULATE_SOLUTION == 1
        testCase.calcSolution(timeManager.getStepSize(),
                              newTimeIdx, advectionManager);
#endif
        timeManager.advance();
        oldTimeIdx.shift();
        advectionManager.output(o1.run("%4.4d", timeManager.getNumStep()), oldTimeIdx);
    }

    return 0;
}
