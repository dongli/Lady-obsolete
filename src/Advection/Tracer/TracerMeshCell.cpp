#include "TracerMeshCell.h"

namespace lady {

TracerMeshCell::TracerMeshCell() {
    x = NULL;
    numConnectedTracer = 0;
    numContainedTracer = 0;
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
#ifndef NDEBUG
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

void TracerMeshCell::disconnect(Tracer *tracer) {
#ifndef NDEBUG
    int i = 0;
    for (; i < numConnectedTracer; ++i) {
        if (connectedTracers[i] == tracer) {
            break;
        }
    }
    assert(i != numConnectedTracer);
#endif
    for (int i = 0; i < numConnectedTracer; ++i) {
        if (connectedTracers[i] == tracer) {
            totalRemapWeight -= remapWeights[i];
            for (int j = i+1; j < numConnectedTracer; ++j) {
                connectedTracers[j-1] = connectedTracers[j];
            }
            numConnectedTracer--;
            break;
        }
    }
}

double TracerMeshCell::getRemapWeight(Tracer *tracer) const {
    for (int i = 0; i < numConnectedTracer; ++i) {
        if (connectedTracers[i] == tracer) {
            return remapWeights[i];
        }
    }
    REPORT_ERROR("Tracer is not connected!");
}

void TracerMeshCell::resetContainedTracers() {
    numContainedTracer = 0;
}

void TracerMeshCell::contain(Tracer *tracer) {
#ifndef NDEBUG
    for (int i = 0; i < numContainedTracer; ++i) {
        if (containedTracers[i] == tracer) {
            REPORT_ERROR("Tracer (ID = " << tracer->getID() <<
                         ") has already been contained!");
        }
    }
#endif
    if (numContainedTracer == containedTracers.size()) {
        containedTracers.push_back(tracer);
    } else {
        containedTracers[numContainedTracer] = tracer;
    }
    tracer->setHostCell(this);
    numContainedTracer++;
}


}