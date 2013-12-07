#include "Tracer.h"

namespace lady {

Tracer::Tracer(int numDim) {
    for (int i = 0; i < q.getNumLevel(); ++i) {
        q.getLevel(i) = new LADY_SPACE_COORD(numDim);
        H.getLevel(i) = new LADY_MATRIX(numDim, numDim);
    }
}

Tracer::~Tracer() {
    for (int i = 0; i < q.getNumLevel(); ++i) {
        delete q.getLevel(i);
        delete H.getLevel(i);
    }
}

Tracer& Tracer::operator=(const Tracer &other) {
    if (this != &other) {
        for (int i = 0; i < q.getNumLevel(); ++i) {
            *(q.getLevel(i)) = *(other.q.getLevel(i));
            *(H.getLevel(i)) = *(other.H.getLevel(i));
        }
    }
    return *this;
}

LADY_SPACE_COORD& Tracer::getX(int timeLevel) {
    return *(q.getLevel(timeLevel));
}

LADY_MATRIX& Tracer::getH(int timeLevel) {
    return *(H.getLevel(timeLevel));
}

}
