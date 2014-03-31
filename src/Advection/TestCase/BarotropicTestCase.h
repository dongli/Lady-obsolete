#ifndef __Lady_BarotropicTestCase__
#define __Lady_BarotropicTestCase__

#include "AdvectionTestCase.h"
#include "barotropic_model.h"

namespace lady {

class BarotropicTestCase : public AdvectionTestCase {
protected:
    barotropic_model::ToyTestCase testCase;
    barotropic_model::BarotropicModel_A_ImplicitMidpoint model;
    geomtk::IOManager<geomtk::RLLDataFile> io;
    int fileIdx;
public:
    BarotropicTestCase();
    virtual ~BarotropicTestCase();

    virtual void init(const geomtk::TimeManager &timeManager);

    Time getStartTime() const;
    Time getEndTime() const;
    double getStepSize() const;

    virtual const LADY_DOMAIN& getDomain() const { return model.getDomain(); }
    virtual const LADY_MESH& getMesh() const { return model.getMesh(); }
    
    void calcInitCond(AdvectionManager &advectionManager);
    void advance(double time, const TimeLevelIndex<2> &timeIdx);
};

}

#endif