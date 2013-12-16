#ifndef __Lady_SolidRotationTestCase__
#define __Lady_SolidRotationTestCase__

#include "AdvectionTestCase.h"

namespace lady {

class SolidRotationTestCase : public AdvectionTestCase {
protected:
    double angleSpeed;
    double U0;
    double alpha;
    double R;
    double H0;
public:
    SolidRotationTestCase();
    virtual ~SolidRotationTestCase();

    Time getStartTime() const;
    Time getEndTime() const;
    double getStepSize() const;

    void advance(double time, int timeLevel);
};

}

#endif