#ifndef __Lady_AdvectionTestCase_h__
#define __Lady_AdvectionTestCase_h__

#include "lady_commons.h"

namespace lady {

class AdvectionTestCase {
protected:
    LADY_DOMAIN *domain;
    LADY_MESH *mesh;
    LADY_VELOCITY_FIELD *v;
public:
    AdvectionTestCase();
    virtual ~AdvectionTestCase();

    const LADY_DOMAIN& getDomain() const;
    const LADY_MESH& getMesh() const;
    const LADY_VELOCITY_FIELD& getVelocityField() const;
    void outputVelocity(const string &fileName, int timeLevel) const;

    /**
     *  Return the start time of the test case.
     *
     *  @return A Time object.
     */
    virtual Time getStartTime() const = 0;

    /**
     *  Return the end time of the test case.
     *
     *  @return A Time object.
     */
    virtual Time getEndTime() const = 0;

    /**
     *  Return the time step size of the test case.
     *
     *  @return The step size in seconds.
     */
    virtual double getStepSize() const = 0;

    /**
     *  Advance the test case one time step.
     *
     *  @param time      the current time in seconds.
     *  @param timeLevel the time level (old or new).
     *
     *  @return None.
     */
    virtual void advance(double time, int timeLevel) = 0;
};

}

#endif
