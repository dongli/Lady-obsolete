#include "Tracer.h"
#include "TracerSkeleton.h"
#include "TracerMeshCell.h"

namespace lady {

Tracer::Tracer(int numDim) : Parcel(numDim) {
    for (int l = 0; l < idx.getNumLevel(); ++l) {
        idx.getLevel(l) = new LADY_MESH_INDEX(numDim);
    }
    skeleton = new TracerSkeleton(this, numDim);
}

Tracer::~Tracer() {
    for (int l = 0; l < idx.getNumLevel(); ++l) {
        delete idx.getLevel(l);
    }
    delete skeleton;
}
    
void Tracer::addSpecies() {
    for (int l = 0; l < ms.getNumLevel(); ++l) {
        ms.getLevel(l).push_back(0.0);
    }
}

double& Tracer::getSpeciesMass(const TimeLevelIndex<2> &timeIdx, int i) {
#ifdef DEBUG
    if (i >= ms.getLevel(timeIdx).size()) {
        REPORT_ERROR("Species index " << i << " exceeds range [0," <<
                     ms.getLevel(timeIdx).size()-1 << "]!");
    }
#endif
    return ms.getLevel(timeIdx)[i];
}

double Tracer::getSpeciesMass(const TimeLevelIndex<2> &timeIdx, int i) const {
#ifdef DEBUG
    if (i >= ms.getLevel(timeIdx).size()) {
        REPORT_ERROR("Species index " << i << " exceeds range [0," <<
                     ms.getLevel(timeIdx).size()-1 << "]!");
    }
#endif
    return ms.getLevel(timeIdx)[i];
}

Tracer& Tracer::operator=(const Tracer &other) {
    Parcel::operator=(other);
    if (this != &other) {
        for (int l = 0; l < idx.getNumLevel(); ++l) {
            *(idx.getLevel(l)) = *(other.idx.getLevel(l));
        }
        *skeleton = *(other.skeleton);
    }
    return *this;
}

LADY_MESH_INDEX& Tracer::getMeshIndex(const TimeLevelIndex<2> &timeIdx) {
    return *(idx.getLevel(timeIdx));
}

TracerSkeleton& Tracer::getSkeleton() {
    return *skeleton;
}

void Tracer::resetConnectedCells() {
    cells.clear();
    totalRemapWeight = 0.0;
}

void Tracer::connect(TracerMeshCell *cell) {
    cells.push_back(cell);
    totalRemapWeight += cell->getWeight(this);
}

const list<TracerMeshCell*>& Tracer::getConnectedCells() const {
    return cells;
}

double Tracer::getTotalRemapWeight() const {
    assert(totalRemapWeight != 0.0);
    return totalRemapWeight;
}

void Tracer::updateDeformationMatrix(const LADY_DOMAIN &domain,
                                     const TimeLevelIndex<2> &timeIdx) {
    // fit "linear" deformation matrix to skeleton
    vector<LADY_BODY_COORD*> &ys = skeleton->getYs();
    vector<LADY_SPACE_COORD*> &xs = skeleton->getXs(timeIdx);
    vec dx0, dx1, dx2, dx3;
    double a1, a2, b1, b2, c1, c2, d1, d2;
    if (idx.getLevel(timeIdx)->isOnPole()) {
        // TODO: Make a tag to indicate whether PS coordinate is set or not.
        q.getLevel(timeIdx)->transformToPS(domain); //>! ensure q is ready
        for (int i = 0; i < xs.size(); ++i) {
            xs[i]->transformToPS(domain);
        }
        dx0 = xs[0]->getPSCoord()-q.getLevel(timeIdx)->getPSCoord();
        dx1 = xs[2]->getPSCoord()-q.getLevel(timeIdx)->getPSCoord();
        dx2 = xs[1]->getPSCoord()-q.getLevel(timeIdx)->getPSCoord();
        dx3 = xs[3]->getPSCoord()-q.getLevel(timeIdx)->getPSCoord();
    } else {
        dx0 = domain.diffCoord(*xs[0], *q.getLevel(timeIdx));
        dx1 = domain.diffCoord(*xs[2], *q.getLevel(timeIdx));
        dx2 = domain.diffCoord(*xs[1], *q.getLevel(timeIdx));
        dx3 = domain.diffCoord(*xs[3], *q.getLevel(timeIdx));
    }
    a1 = dx0(0)/(*ys[0])(0);
    a2 = dx1(0)/(*ys[2])(0);
    b1 = dx2(0)/(*ys[1])(1);
    b2 = dx3(0)/(*ys[3])(1);
    c1 = dx0(1)/(*ys[0])(0);
    c2 = dx1(1)/(*ys[2])(0);
    d1 = dx2(1)/(*ys[1])(1);
    d2 = dx3(1)/(*ys[3])(1);
    (*H.getLevel(timeIdx))(0, 0) = (a1+a2)*0.5;
    (*H.getLevel(timeIdx))(0, 1) = (b1+b2)*0.5;
    (*H.getLevel(timeIdx))(1, 0) = (c1+c2)*0.5;
    (*H.getLevel(timeIdx))(1, 1) = (d1+d2)*0.5;
    // update inversion and determinant
    detH.getLevel(timeIdx) = det(*H.getLevel(timeIdx));
    assert(detH.getLevel(timeIdx) > 0.0);
    *invH.getLevel(timeIdx) = inv(*H.getLevel(timeIdx));
    updateShapeSize(domain, timeIdx);
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
