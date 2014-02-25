#ifndef __Lady_DeformMatrixFitting__
#define __Lady_DeformMatrixFitting__

#include "lady_commons.h"
#include "Tracer.h"

namespace lady {

/**
 *  This struct includes the necessary information for optimization. Note the
 *  center is shifted to origin, so x should be processed accordingly.
 */
struct Data {
    const LADY_DOMAIN *domain;
    vector<vec> x, y;
    int count;
    bool debug;
};

/**
 *  This class is used to fit the deformation matrix from tracer skeleton.
 */
class DeformMatrixFitting {
    Data data;
    nlopt::opt *a;
public:
    DeformMatrixFitting(const LADY_DOMAIN &domain, Tracer *tracer);
    ~DeformMatrixFitting();
    void fit(const TimeLevelIndex<2> &timeIdx, Tracer *tracer);
#ifdef DEBUG
    static double maxObjectiveValue;
#endif
private:
    /**
     *  Objective function for optimization.
     *
     *  @param n     the parameter number.
     *  @param x     the parameters.
     *  @param grad  the gradient of function on parameters.
     *  @param data_ the extra data.
     *
     *  @return The objective function value.
     */
    static double objective(unsigned int n,  const double *x, double *grad,
                            void *data_);
    
    /**
     *  Nonlinear inequality constraint for optimization.
     *
     *  @param n     the parameter number.
     *  @param x     the parameters.
     *  @param grad  the gradient of function on parameters.
     *  @param data_ the extra data.
     *
     *  @return The constraint value.
     */
    static double constraint(unsigned int n,  const double *x, double *grad,
                             void *data_);
};

}

#endif
