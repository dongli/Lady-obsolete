#ifndef __Lady_SolidRotationTestCase__
#define __Lady_SolidRotationTestCase__

#include "AdvectionTestCase.h"

namespace lady {

class SolidRotationTestCase : public AdvectionTestCase {
protected:
    double angleSpeed, U0, alpha;
    LADY_SPACE_COORD *axisPole, *c0, *cr0;
    double R, H0;
public:
    SolidRotationTestCase();
    virtual ~SolidRotationTestCase();

    Time getStartTime() const;
    Time getEndTime() const;
    double getStepSize() const;

    void calcInitCond(AdvectionManager &advectionManager);
    void calcSolution(double time, LADY_SCALAR_FIELD &q);
    void advance(double time, int timeLevel);
};

}

#endif