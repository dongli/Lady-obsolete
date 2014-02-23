#include "TracerMeshCell.h"

namespace lady {

TracerMeshCell::TracerMeshCell() {
    x = NULL;
}

TracerMeshCell::~TracerMeshCell() {
    if (x != NULL) {
        delete x;
    }
}

void TracerMeshCell::resetSpeciesMass() {
    for (int i = 0; i < m.size(); ++i) {
        m[i] = 0.0;
    }
}

void TracerMeshCell::resetConnectedTracers() {
    tracers.clear();
    totalRemapWeight = 0.0;
}

void TracerMeshCell::connect(Tracer *tracer, double weight) {
#ifdef DEBUG
    if (tracers.count(tracer) != 0) {
        REPORT_ERROR("Tracer (ID = " << tracer->getID() <<
                     ") has already been connected!");
    }
#endif
    tracers[tracer] = weight;
    totalRemapWeight += weight;
}

double TracerMeshCell::getWeight(Tracer *tracer) const {
#ifdef DEBUG
    if (tracers.count(tracer) == 0) {
        REPORT_ERROR("Tracer (ID = " << tracer->getID() <<
                     ") has not yet been connected!");
    }
#endif
    return tracers.at(tracer);
}

}