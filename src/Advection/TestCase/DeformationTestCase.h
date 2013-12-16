#ifndef __Lady_DeformationTestCase__
#define __Lady_DeformationTestCase__

#include "AdvectionTestCase.h"

namespace lady {

class DeformationTestCase : public AdvectionTestCase {
public:
    enum SubCase {
        CASE1, CASE2, CASE3, CASE4
    };
    enum InitCond {
        COSINE_HILL, GAUSSIAN_HILL, SLOTTED_CYLINDERS
    };
protected:
    double T;
    SubCase subCase;
    InitCond initCond;
public:
    DeformationTestCase(SubCase subCase, InitCond initCond);
    ~DeformationTestCase();

    Time getStartTime() const;
    Time getEndTime() const;
    double getStepSize() const;

    void advance(double time, int timeLevel);
};

}

#endif