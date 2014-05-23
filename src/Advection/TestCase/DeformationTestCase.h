#ifndef __Lady_DeformationTestCase__
#define __Lady_DeformationTestCase__

#include "AdvectionTestCase.h"

namespace lady {

class DeformationTestCase : public AdvectionTestCase {
public:
    enum SubCase {
        CASE1, CASE2, CASE3, CASE4
    };
protected:
    double period;
    SubCase subCase;
public:
    DeformationTestCase(SubCase subCase);
    ~DeformationTestCase();

    virtual void init(const geomtk::TimeManager &timeManager);

    Time getStartTime() const;
    Time getEndTime() const;
    double getStepSize() const;

    void calcInitCond(AdvectionManager &advectionManager);
    void advance(double time, const TimeLevelIndex<2> &timeIdx);
};

}

#endif