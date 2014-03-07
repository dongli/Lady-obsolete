#include "Parcel.h"
#include "ShapeFunction.h"

namespace lady {

Parcel::Parcel(int numDim) {
    for (int l = 0; l < q.getNumLevel(); ++l) {
        q.getLevel(l) = new LADY_SPACE_COORD(numDim);
        H.getLevel(l) = new LADY_MATRIX(numDim, numDim);
        invH.getLevel(l) = new LADY_MATRIX(numDim, numDim);
    }
}

Parcel::~Parcel() {
    for (int l = 0; l < q.getNumLevel(); ++l) {
        delete q.getLevel(l);
        delete H.getLevel(l);
        delete invH.getLevel(l);
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
        }
    }
    return *this;
}

int Parcel::getID() const {
    return ID;
}

void Parcel::setID(int ID) {
    this->ID = ID;
}

LADY_SPACE_COORD& Parcel::getX(const TimeLevelIndex<2> &timeIdx) {
    return *(q.getLevel(timeIdx));
}

LADY_MATRIX& Parcel::getH(const TimeLevelIndex<2> &timeIdx) {
    return *(H.getLevel(timeIdx));
}
    
LADY_MATRIX& Parcel::getInvH(const TimeLevelIndex<2> &timeIdx) {
    return *(invH.getLevel(timeIdx));
}
    
void Parcel::getSpaceCoord(const LADY_DOMAIN &domain,
                           const TimeLevelIndex<2> &timeIdx,
                           const LADY_BODY_COORD &y, LADY_SPACE_COORD &x) {
    if (LADY_IS_SPHERE_DOMAIN) {
        // In sphere domain, we calculate deformation matrix stuffs on local
        // stereographic projection of tracer centroid.
        x() = (*H.getLevel(timeIdx))*y();
        domain.projectBack(geomtk::SphereDomain::STEREOGRAPHIC,
                           *q.getLevel(timeIdx), x, x());
    } else {
        REPORT_ERROR("Under construction!");
        x() = (*q.getLevel(timeIdx))()+(*H.getLevel(timeIdx))*y();
    }
}

void Parcel::getBodyCoord(const LADY_DOMAIN &domain,
                          const TimeLevelIndex<2> &timeIdx,
                          const LADY_SPACE_COORD &x, LADY_BODY_COORD &y) {
    if (LADY_IS_SPHERE_DOMAIN) {
        domain.project(geomtk::SphereDomain::STEREOGRAPHIC,
                       *q.getLevel(timeIdx), x, y());
        y() = (*invH.getLevel(timeIdx))*y();
    } else {
        REPORT_ERROR("Under construction!");
        y() = (*invH.getLevel(timeIdx))*(x()-(*q.getLevel(timeIdx))());
    }
}
    
double Parcel::getShapeFunction(const TimeLevelIndex<2> &timeIdx,
                                const LADY_BODY_COORD &y) {
    double f;
    ShapeFunction::evalFunc(y, f);
    f /= detH.getLevel(timeIdx);
    return f;
}

void Parcel::updateShapeSize(const LADY_DOMAIN &domain,
                             const TimeLevelIndex<2> &timeIdx) {
    LADY_BODY_COORD y(domain.getNumDim());
    LADY_SPACE_COORD x(domain.getNumDim());
    if (domain.getNumDim() == 2) {
        y(0) = 1.0; y(1) = 0.0;
        getSpaceCoord(domain, timeIdx, y, x);
        double d1 = domain.calcDistance(x, *(q.getLevel(timeIdx)));
        y(0) = 0.0; y(1) = 1.0;
        getSpaceCoord(domain, timeIdx, y, x);
        double d2 = domain.calcDistance(x, *(q.getLevel(timeIdx)));
        shapeSize.getLevel(timeIdx)(0) = fmax(d1, d2);
        shapeSize.getLevel(timeIdx)(1) = fmin(d1, d2);
    } else {
        REPORT_ERROR("Under construction!");
    }
}

const vec::fixed<2>& Parcel::getShapeSize(const TimeLevelIndex<2> &timeIdx) const {
    return shapeSize.getLevel(timeIdx);
}

}
