#include "Tracer.h"

namespace lady {

Tracer::Tracer(int numDim) {
    for (int i = 0; i < q.getNumLevel(); ++i) {
        q.getLevel(i) = new LADY_SPACE_COORD(numDim);
        H.getLevel(i) = new LADY_MATRIX(numDim, numDim);
        idx.getLevel(i) = new LADY_MESH_INDEX(numDim);
    }
}

Tracer::~Tracer() {
    for (int i = 0; i < q.getNumLevel(); ++i) {
        delete q.getLevel(i);
        delete H.getLevel(i);
        delete idx.getLevel(i);
    }
}

Tracer& Tracer::operator=(const Tracer &other) {
    if (this != &other) {
        for (int i = 0; i < q.getNumLevel(); ++i) {
            *(q.getLevel(i)) = *(other.q.getLevel(i));
            *(H.getLevel(i)) = *(other.H.getLevel(i));
            *(idx.getLevel(i)) = *(other.idx.getLevel(i));
        }
    }
    return *this;
}

int Tracer::getID() const {
    return ID;
}

void Tracer::setID(int ID) {
    this->ID = ID;
}

LADY_SPACE_COORD& Tracer::getX(int timeLevel) {
    return *(q.getLevel(timeLevel));
}

LADY_MATRIX& Tracer::getH(int timeLevel) {
    return *(H.getLevel(timeLevel));
}

LADY_MESH_INDEX& Tracer::getMeshIndex(int timeLevel) {
    return *(idx.getLevel(timeLevel));
}

}
