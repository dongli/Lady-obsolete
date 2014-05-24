#include "lady.h"

int main(int argc, const char *argv[])
{
    geomtk::ConfigManager configManager;
    lady::AdvectionTestCase *testCase;
    lady::AdvectionManager advectionManager;
    geomtk::TimeManager timeManager;
    geomtk::TimeLevelIndex<2> oldTimeIdx;
    geomtk::StampString *o1;
    // -------------------------------------------------------------------------
    // parse configuration
    configManager.parse(argv[1]);
    // -------------------------------------------------------------------------
    // choose test case
    bool isTrueSolution = false;
    std::string testCaseName, subcaseName = "";
    configManager.getValue("lasm", "test_case", testCaseName);
    if (testCaseName == "rotation") {
        testCase = new lady::SolidRotationTestCase();
        if (configManager.hasKey("lasm", "is_true_solution")) {
            configManager.getValue("lasm", "is_true_solution", isTrueSolution);
        }
    } else if (testCaseName == "deform") {
        testCase = new lady::DeformationTestCase();
        if (configManager.hasKey("lasm", "sub_case")) {
            configManager.getValue("lasm", "sub_case", subcaseName);
            testCase->selectSubcase(subcaseName);
        }
    } else if (testCaseName == "barotropic") {
        testCase = new lady::BarotropicTestCase();
    } else {
        REPORT_ERROR("Unknown test_case \"" << testCaseName << "\"!");
    }
    std::string outputPrefix = "";
    if (configManager.hasKey("lasm", "output_prefix")) {
        configManager.getValue("lasm", "output_prefix", outputPrefix);
    } else {
        outputPrefix = "tracers."+testCaseName+".";
        if (isTrueSolution) {
            outputPrefix += "true.";
        }
        int nx, ny;
        configManager.getValue("lasm", "num_parcel_x", nx);
        configManager.getValue("lasm", "num_parcel_y", ny);
        if (subcaseName == "") {
            outputPrefix += std::to_string(nx)+"x"+std::to_string(ny);
        } else {
            outputPrefix += subcaseName+"."+std::to_string(nx)+"x"+std::to_string(ny);
        }
    }
    o1 = new geomtk::StampString(outputPrefix+".", ".nc");
    // -------------------------------------------------------------------------
    // initialization
    timeManager.init(testCase->getStartTime(), testCase->getEndTime(),
                     testCase->getStepSize());
    testCase->init(timeManager);
    advectionManager.init(testCase->getDomain(), testCase->getMesh(), configManager);
    
    testCase->calcInitCond(advectionManager);
    advectionManager.output(o1->run("%4.4d", timeManager.getNumStep()), oldTimeIdx);
    testCase->advance(timeManager.getSeconds(), oldTimeIdx);
    // -------------------------------------------------------------------------
    // integration loop
    while (!timeManager.isFinished()) {
        geomtk::TimeLevelIndex<2> newTimeIdx = oldTimeIdx+1;
        double time = timeManager.getSeconds()+timeManager.getStepSize();
        testCase->advance(time, newTimeIdx);
        advectionManager.advance(timeManager.getStepSize(), newTimeIdx,
                                 testCase->getVelocityField());
        if (isTrueSolution) {
            testCase->calcSolution(timeManager.getStepSize(),
                                   newTimeIdx, advectionManager);
        }
        timeManager.advance();
        oldTimeIdx.shift();
        advectionManager.output(o1->run("%4.4d", timeManager.getNumStep()), oldTimeIdx);
    }
    delete testCase;
    delete o1;
    return 0;
}
