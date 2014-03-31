#include "TracerMeshCell.h"

namespace lady {

TracerMeshCell::TracerMeshCell() {
    x = NULL;
    numConnectedTracer = 0;
}

TracerMeshCell::~TracerMeshCell() {
    if (x != NULL) {
        delete x;
    }
}

void TracerMeshCell::resetConnectedTracers() {
    numConnectedTracer = 0;
    totalRemapWeight = 0;
}

void TracerMeshCell::connect(Tracer *tracer, double weight) {
#ifdef DEBUG
    for (int i = 0; i < numConnectedTracer; ++i) {
        if (connectedTracers[i] == tracer) {
            REPORT_ERROR("Tracer (ID = " << tracer->getID() <<
                         ") has already been connected!");
        }
    }
    if (numConnectedTracer == 0) {
        assert(totalRemapWeight == 0);
    }
#endif
    if (numConnectedTracer == connectedTracers.size()) {
        connectedTracers.push_back(tracer);
        remapWeights.push_back(weight);
    } else {
        connectedTracers[numConnectedTracer] = tracer;
        remapWeights[numConnectedTracer] = weight;
    }
    numConnectedTracer++;
    totalRemapWeight += weight;
}

double TracerMeshCell::getRemapWeight(Tracer *tracer) const {
    for (int i = 0; i < numConnectedTracer; ++i) {
        if (connectedTracers[i] == tracer) {
            return remapWeights[i];
        }
    }
    REPORT_ERROR("Tracer is not connected!");
}

}