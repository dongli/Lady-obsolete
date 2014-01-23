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
    // TODO: How to hide sphere domain details?
    if (idx.getLevel(timeIdx)->isOnPole()) {
        x.getPSCoord() = q.getLevel(timeIdx)->getPSCoord()+(*H.getLevel(timeIdx))*y();
        x.transformFromPS(domain, idx.getLevel(timeIdx)->getPole());
    } else {
        x() = (*q.getLevel(timeIdx))()+(*H.getLevel(timeIdx))*y();
        domain.check(x);
    }
}

void Parcel::getBodyCoord(const LADY_DOMAIN &domain,
                          const TimeLevelIndex<2> &timeIdx,
                          const LADY_SPACE_COORD &x, LADY_BODY_COORD &y) {
    // TODO: How to hide sphere domain details?
    if (idx.getLevel(timeIdx)->isOnPole()) {
        y() = (*invH.getLevel(timeIdx))*(x.getPSCoord()-
                                         q.getLevel(timeIdx)->getPSCoord());
    } else {
        vec dx = domain.diffCoord(x, *q.getLevel(timeIdx));
        y() = (*invH.getLevel(timeIdx))*dx;
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
    y(0) = 1.0;
    y(1) = 0.0;
    getSpaceCoord(domain, timeIdx, y, x);
    double d1 = domain.calcDistance(x, *(q.getLevel(timeIdx)));
    y(0) = 0.0;
    y(1) = 1.0;
    getSpaceCoord(domain, timeIdx, y, x);
    double d2 = domain.calcDistance(x, *(q.getLevel(timeIdx)));
    shapeSize.getLevel(timeIdx)(0) = fmax(d1, d2);
    shapeSize.getLevel(timeIdx)(1) = fmin(d1, d2);
}

const vec::fixed<2>& Parcel::getShapeSize(const TimeLevelIndex<2> &timeIdx) const {
    return shapeSize.getLevel(timeIdx);
}

}