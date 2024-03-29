#ifndef __Lady_AdvectionTestCase_h__
#define __Lady_AdvectionTestCase_h__

#include "lady_commons.h"
#include "AdvectionManager.h"

namespace lady {

class AdvectionTestCase {
protected:
    LADY_DOMAIN *domain;
    LADY_MESH *mesh;
    const geomtk::TimeManager *timeManager;
    LADY_VELOCITY_FIELD velocity;
    vector<LADY_SCALAR_FIELD*> q;
public:
    AdvectionTestCase();
    virtual ~AdvectionTestCase();

    virtual void init(const geomtk::TimeManager &timeManager);

    virtual const LADY_DOMAIN& getDomain() const { return *domain; }
    virtual const LADY_MESH& getMesh() const { return *mesh; }
    virtual const LADY_VELOCITY_FIELD& getVelocityField() const { return velocity; }
    
    /**
     *  Output velocity field.
     *
     *  @param fileName   the output file name.
     *  @param oldTimeIdx the old time level index.
     */
    void outputVelocity(const string &fileName,
                        const TimeLevelIndex<2> &oldTimeIdx) const;

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
     *  Calculate initial condition and set tracers. This base function just
     *  transfer the settings of derived ones to tracers, so it will be called
     *  at the end of derived ones.
     */
    virtual void calcInitCond(AdvectionManager &advectionManager);

    /**
     *  Advance the test case one time step.
     *
     *  @param time    the time in seconds.
     *  @param timeIdx the time level index.
     */
    virtual void advance(double time, const TimeLevelIndex<2> &timeIdx) = 0;

    /**
     *  Calculate the solution of the test case if any, and reset the tracers
     *  for latter outputting.
     *
     *  @param time             the time in seconds.
     *  @param timeIdx          the time level index.
     *  @param advectionManager the advection manager.
     */
    virtual void calcSolution(double time, const TimeLevelIndex<2> &timeIdx,
                              AdvectionManager &advectionManager);
protected:
    /**
     *  Calculate the solution of the test case if any.
     *
     *  @param time    the time in seconds.
     *  @param timeIdx the time level index.
     *  @param q       the output solution.
     */
    virtual void calcSolution(double time, const TimeLevelIndex<2> &timeIdx,
                              LADY_SCALAR_FIELD &q);
};

}

#endif
