#include <TracerManager.h>

namespace lady
{

TracerManager::TracerManager() {
}

TracerManager::~TracerManager() {
}

void TracerManager::init(const LADY_DOMAIN &domain, int numTracer) {
    tracers.resize(numTracer);
    for (LADY_LIST<Tracer*>::iterator t = tracers.begin(); t != tracers.end(); ++t) {
        *t = new Tracer(domain.getNumDim());
    }
    REPORT_NOTICE(numTracer << " tracers are initialized.");
}

}
