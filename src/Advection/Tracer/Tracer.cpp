#include "Tracer.h"
#include "TracerSkeleton.h"
#include "TracerMeshCell.h"
#include "DeformMatrixFitting.h"

namespace lady {

Tracer::Tracer(int numDim) : Parcel(numDim) {
    for (int l = 0; l < idx.getNumLevel(); ++l) {
        idx.getLevel(l) = new LADY_MESH_INDEX(numDim);
    }
    skeleton = new TracerSkeleton(this, numDim);
    deformMatrixFitting = NULL;
}

Tracer::~Tracer() {
    for (int l = 0; l < idx.getNumLevel(); ++l) {
        delete idx.getLevel(l);
    }
    delete skeleton;
    if (deformMatrixFitting != NULL) {
        delete deformMatrixFitting;
    }
}
    
void Tracer::addSpecies() {
    m.push_back(0.0);
}

double& Tracer::getSpeciesMass(int s) {
#ifdef DEBUG
    if (s >= m.size()) {
        REPORT_ERROR("Species index " << s << " exceeds range [0," <<
                     m.size()-1 << "]!");
    }
#endif
    return m[s];
}

double Tracer::getSpeciesMass(int s) const {
#ifdef DEBUG
    if (s >= m.size()) {
        REPORT_ERROR("Species index " << s << " exceeds range [0," <<
                     m.size()-1 << "]!");
    }
#endif
    return m[s];
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

void Tracer::resetSpeciesMass() {
    for (int s = 0; s < m.size(); ++s) {
        m[s] = 0.0;
    }
}

void Tracer::connect(TracerMeshCell *cell, double weight) {
    cells.push_back(cell);
    totalRemapWeight += weight;
}

list<TracerMeshCell*>& Tracer::getConnectedCells() {
    return cells;
}

double Tracer::getTotalRemapWeight() const {
    assert(totalRemapWeight != 0.0);
    return totalRemapWeight;
}

void Tracer::updateDeformMatrix(const LADY_DOMAIN &domain,
                                const TimeLevelIndex<2> &timeIdx,
                                bool isFirstTime) {
    if (isFirstTime) {
        deformMatrixFitting = new DeformMatrixFitting(domain, this);
        *q.getLevel(1) = *q.getLevel(0);
        *idx.getLevel(1) = *idx.getLevel(0);
        for (int i = 0; i < skeleton->getYs().size(); ++i) {
            *skeleton->getXs(timeIdx+1)[i] = *skeleton->getXs(timeIdx)[i];
        }
        deformMatrixFitting->fit(timeIdx+1, this);
        *H.getLevel(0) = *H.getLevel(1);
    } else {
        deformMatrixFitting->fit(timeIdx, this);
    }
    detH.getLevel(timeIdx) = det(*H.getLevel(timeIdx));
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
