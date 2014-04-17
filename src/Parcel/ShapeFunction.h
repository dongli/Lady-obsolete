#ifndef __Lady_ShapeFunction__
#define __Lady_ShapeFunction__

#include "lady_commons.h"

namespace lady {

class ShapeFunction {
public:
    static double J;
    // one-dimensional quadrature nodes and weights
    static vec nodes;
    static vec weights;
private:
    static const LADY_DOMAIN *domain;
    static double maxValue;
public:
    static void init(const LADY_DOMAIN &domain);
    static double getMaxValue() { return maxValue; }
    static void evalFunc(const LADY_BODY_COORD& y, double &f);
    static void evalDerv(const LADY_BODY_COORD& y, vec &d);

private:
    ShapeFunction() {}
    virtual ~ShapeFunction() {}
};

}

#endif