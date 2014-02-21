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

void TracerMeshCell::setCoord(const LADY_SPACE_COORD &x) {
    this->x = new LADY_SPACE_COORD(x);
}

const LADY_SPACE_COORD& TracerMeshCell::getCoord() const {
    return *x;
}

void TracerMeshCell::setVolume(double volume) {
    this->volume = volume;
}

double TracerMeshCell::getVolume() const {
    return volume;
}

void TracerMeshCell::addSpecies() {
    m.push_back(0.0);
}

double& TracerMeshCell::getSpeciesMass(int speciesIdx) {
    return m[speciesIdx];
}

double TracerMeshCell::getSpeciesMass(int speciesIdx) const {
    return m[speciesIdx];
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

double TracerMeshCell::getTotalRemapWeight() const {
    assert(totalRemapWeight != 0.0);
    return totalRemapWeight;
}

}