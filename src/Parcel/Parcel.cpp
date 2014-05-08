#include "Parcel.h"
#include "ShapeFunction.h"

namespace lady {

Parcel::Parcel(int numDim) {
    for (int l = 0; l < q.getNumLevel(); ++l) {
        q.getLevel(l) = new LADY_SPACE_COORD(numDim);
        H.getLevel(l) = new LADY_MATRIX(numDim, numDim);
        invH.getLevel(l) = new LADY_MATRIX(numDim, numDim);
        idx.getLevel(l) = new LADY_MESH_INDEX(numDim);
    }
}

Parcel::~Parcel() {
    for (int l = 0; l < q.getNumLevel(); ++l) {
        delete q.getLevel(l);
        delete H.getLevel(l);
        delete invH.getLevel(l);
        delete idx.getLevel(l);
    }
}
    
Parcel& Parcel::operator=(const Parcel &other) {
    if (this != &other) {
        for (int l = 0; l < q.getNumLevel(); ++l) {
            *(q.getLevel(l)) = *(other.q.getLevel(l));
            *(H.getLevel(l)) = *(other.H.getLevel(l));
            detH.getLevel(l) = other.detH.getLevel(l);
            *(invH.getLevel(l)) = *(other.invH.getLevel(l));
            shapeSize.getLevel(l) = other.shapeSize.getLevel(l);
            *(idx.getLevel(l)) = *(other.idx.getLevel(l));
        }
    }
    return *this;
}
    
void Parcel::getSpaceCoord(const LADY_DOMAIN &domain,
                           const TimeLevelIndex<2> &timeIdx,
                           const LADY_BODY_COORD &y, LADY_SPACE_COORD &x) {
#ifdef LADY_USE_SPHERE_DOMAIN
    // In sphere domain, we calculate deformation matrix stuffs on local
    // stereographic projection of tracer centroid.
    x() = (*H.getLevel(timeIdx))*y();
    domain.projectBack(geomtk::SphereDomain::STEREOGRAPHIC,
                       *q.getLevel(timeIdx), x, x());
#else
    REPORT_ERROR("Under construction!");
    x() = (*q.getLevel(timeIdx))()+(*H.getLevel(timeIdx))*y();
#endif
}

void Parcel::getBodyCoord(const LADY_DOMAIN &domain,
                          const TimeLevelIndex<2> &timeIdx,
                          const LADY_SPACE_COORD &x, LADY_BODY_COORD &y) {
#ifdef LADY_USE_SPHERE_DOMAIN
    domain.project(geomtk::SphereDomain::STEREOGRAPHIC,
                   *q.getLevel(timeIdx), x, y());
    y() = (*invH.getLevel(timeIdx))*y();
#else
    REPORT_ERROR("Under construction!");
    y() = (*invH.getLevel(timeIdx))*(x()-(*q.getLevel(timeIdx))());
#endif
}
    
double Parcel::getShapeFunction(const TimeLevelIndex<2> &timeIdx,
                                const LADY_BODY_COORD &y) {
    double f;
    ShapeFunction::evalFunc(y, f);
//    f /= detH.getLevel(timeIdx);
    return f;
}

void Parcel::updateShapeSize(const LADY_DOMAIN &domain,
                             const TimeLevelIndex<2> &timeIdx) {
    LADY_BODY_COORD y(domain.getNumDim());
    LADY_SPACE_COORD x(domain.getNumDim());
    for (int m = 0; m < domain.getNumDim(); ++m) {
        y() = (*invH.getLevel(timeIdx))*(*H.getLevel(timeIdx))*V.col(m);
        getSpaceCoord(domain, timeIdx, y, x);
        shapeSize.getLevel(timeIdx)[m] = domain.calcDistance(x, *q.getLevel(timeIdx));
    }
}

}
