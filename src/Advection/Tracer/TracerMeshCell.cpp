#include "TracerMeshCell.h"

namespace lady {

TracerMeshCell::TracerMeshCell() {
    x = NULL;
    maxNumConnectedTracer = 1000;
    numConnectedTracer = 0;
    connectedTracers = new Tracer*[maxNumConnectedTracer];
    remapWeights = new double[maxNumConnectedTracer];
}

TracerMeshCell::~TracerMeshCell() {
    if (x != NULL) {
        delete x;
    }
    delete [] connectedTracers;
    delete [] remapWeights;
}

void TracerMeshCell::resetSpeciesMass() {
    for (int i = 0; i < m.size(); ++i) {
        m[i] = 0.0;
    }
}

void TracerMeshCell::resetConnectedTracers() {
    numConnectedTracer = 0;
    totalRemapWeight = 0;
}

void TracerMeshCell::connect(Tracer *tracer, double weight) {
    if (maxNumConnectedTracer <= numConnectedTracer) {
        REPORT_ERROR("Limit of tracer number (" << maxNumConnectedTracer <<
                     ") has been reached!");
    }
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
    connectedTracers[numConnectedTracer] = tracer;
    remapWeights[numConnectedTracer] = weight;
    numConnectedTracer++;
    totalRemapWeight += weight;
}

double TracerMeshCell::getRemapWeight(Tracer *tracer) const {
    int i;
    for (i = 0; i < numConnectedTracer; ++i) {
        if (connectedTracers[i] == tracer) {
            return remapWeights[i];
        }
    }
    REPORT_ERROR("Tracer is not connected!");
}

}