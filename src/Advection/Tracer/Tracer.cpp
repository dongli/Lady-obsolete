#include "Tracer.h"
#include "TracerSkeleton.h"
#include "ShapeFunction.h"

namespace lady {

Tracer::Tracer(int numDim) {
    for (int l = 0; l < q.getNumLevel(); ++l) {
        q.getLevel(l) = new LADY_SPACE_COORD(numDim);
        H.getLevel(l) = new LADY_MATRIX(numDim, numDim);
        idx.getLevel(l) = new LADY_MESH_INDEX(numDim);
        invH.getLevel(l) = new LADY_MATRIX(numDim, numDim);
    }
    skeleton = new TracerSkeleton(this, numDim);
}

Tracer::~Tracer() {
    for (int l = 0; l < q.getNumLevel(); ++l) {
        delete q.getLevel(l);
        delete H.getLevel(l);
        delete idx.getLevel(l);
        delete invH.getLevel(l);
    }
    delete skeleton;
}
    
void Tracer::addSpecies() {
    for (int l = 0; l < ms.getNumLevel(); ++l) {
        ms.getLevel(l).resize(ms.getLevel(l).size()+1);
    }
}

Tracer& Tracer::operator=(const Tracer &other) {
    if (this != &other) {
        for (int l = 0; l < q.getNumLevel(); ++l) {
            *(q.getLevel(l)) = *(other.q.getLevel(l));
            *(H.getLevel(l)) = *(other.H.getLevel(l));
            *(idx.getLevel(l)) = *(other.idx.getLevel(l));
            detH.getLevel(l) = other.detH.getLevel(l);
            *(invH.getLevel(l)) = *(other.invH.getLevel(l));
        }
        *skeleton = *(other.skeleton);
    }
    return *this;
}

int Tracer::getID() const {
    return ID;
}

void Tracer::setID(int ID) {
    this->ID = ID;
}

LADY_SPACE_COORD& Tracer::getX(const TimeLevelIndex<2> &timeIdx) {
    return *(q.getLevel(timeIdx));
}

LADY_MATRIX& Tracer::getH(const TimeLevelIndex<2> &timeIdx) {
    return *(H.getLevel(timeIdx));
}
    
LADY_MATRIX& Tracer::getInvH(const TimeLevelIndex<2> &timeIdx) {
    return *(invH.getLevel(timeIdx));
}

LADY_MESH_INDEX& Tracer::getMeshIndex(const TimeLevelIndex<2> &timeIdx) {
    return *(idx.getLevel(timeIdx));
}

TracerSkeleton& Tracer::getSkeleton() {
    return *skeleton;
}
    
void Tracer::getSpaceCoord(const LADY_DOMAIN &domain,
                           const TimeLevelIndex<2> &timeIdx,
                           const LADY_BODY_COORD &y, LADY_SPACE_COORD &x) {
    x() = (*q.getLevel(timeIdx))()+(*H.getLevel(timeIdx))*y();
}

void Tracer::getBodyCoord(const LADY_DOMAIN &domain,
                          const TimeLevelIndex<2> &timeIdx,
                          const LADY_SPACE_COORD &x, LADY_BODY_COORD &y) {
    vec dx = domain.diffCoord(x, *q.getLevel(timeIdx));
    y() = (*invH.getLevel(timeIdx))*dx;
}
    
double Tracer::getShapeFunction(const TimeLevelIndex<2> &timeIdx,
                                const LADY_BODY_COORD &y) {
    double f;
    ShapeFunction::evalFunc(y, f);
    f /= detH.getLevel(timeIdx);
    return f;
}
    
void Tracer::updateH(const LADY_DOMAIN &domain,
                     const TimeLevelIndex<2> &timeIdx) {
    // fit "linear" deformation matrix to skeleton
    vector<LADY_BODY_COORD*> &ys = skeleton->getYs();
    vector<LADY_SPACE_COORD*> &xs = skeleton->getXs(timeIdx);
    vec dx0 = domain.diffCoord(*xs[0], *q.getLevel(timeIdx));
    vec dx1 = domain.diffCoord(*xs[1], *q.getLevel(timeIdx));
    vec dx2 = domain.diffCoord(*xs[2], *q.getLevel(timeIdx));
    vec dx3 = domain.diffCoord(*xs[3], *q.getLevel(timeIdx));
    double a1 = dx0(0)/(*ys[0])(0);
    double a2 = dx1(0)/(*ys[1])(0);
    double b1 = dx2(0)/(*ys[2])(1);
    double b2 = dx3(0)/(*ys[3])(1);
    double c1 = dx0(1)/(*ys[0])(0);
    double c2 = dx1(1)/(*ys[1])(0);
    double d1 = dx2(1)/(*ys[2])(1);
    double d2 = dx3(1)/(*ys[3])(1);
    (*H.getLevel(timeIdx))(0, 0) = (a1+a2)*0.5;
    (*H.getLevel(timeIdx))(0, 1) = (b1+b2)*0.5;
    (*H.getLevel(timeIdx))(1, 0) = (c1+c2)*0.5;
    (*H.getLevel(timeIdx))(1, 1) = (d1+d2)*0.5;
    // update inversion and determinant
    detH.getLevel(timeIdx) = det(*H.getLevel(timeIdx));
    *invH.getLevel(timeIdx) = inv(*H.getLevel(timeIdx));
}
    
void Tracer::selfInspect(const LADY_DOMAIN &domain,
                         const TimeLevelIndex<2> &timeIdx) {
    static const double offsetThreshold = 1.0*RAD/domain.getRadius();
    vector<LADY_BODY_COORD*> &ys = skeleton->getYs();
    vector<LADY_SPACE_COORD*> &xs = skeleton->getXs(timeIdx);
    for (int i = 0; i < ys.size(); ++i) {
        LADY_SPACE_COORD x(domain.getNumDim());
        getSpaceCoord(domain, timeIdx, *ys[i], x);
        double offset = domain.calcDistance(*xs[i], x);
        if (offset >= offsetThreshold) {
            cout << ID << endl;
        }
    }
}

}
