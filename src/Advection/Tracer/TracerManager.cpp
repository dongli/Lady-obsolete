#include <TracerManager.h>

namespace lady
{

TracerManager::TracerManager() {
    REPORT_NOTICE("Tracer manager is online");
}

    TracerManager::~TracerManager() {
        REPORT_NOTICE("Tracer manager is offline");
}

void TracerManager::init(const LADY_DOMAIN &domain, const LADY_MESH &mesh,
                         int numTracer) {
    this->domain = &domain;
    tracers.resize(numTracer);
    int ID = 0;
    LADY_LIST<Tracer*>::iterator t = tracers.begin();
    for (; t != tracers.end(); ++t) {
        *t = new Tracer(domain.getNumDim());
        (*t)->setID(ID++);
    }
    // -------------------------------------------------------------------------
    // set initial centroid coordinates and deformation matrices of tracers
    double dx[domain.getNumDim()];
    int nt[domain.getNumDim()];
#define DEBUG_CROSS_POLE 0
#if DEBUG_CROSS_POLE == 1
    assert(dynamic_cast<const geomtk::SphereDomain*>(&domain));
    nt[0] = numTracer; nt[1] = 1;
    dx[0] = PI2/numTracer; dx[1] = 5*RAD;
#else
    int tmp1 = numTracer;
    for (int m = 0; m < domain.getNumDim(); ++m) {
        int tmp2 = ceil(pow(tmp1, 1.0/(domain.getNumDim()-m)));
        for (; tmp2 >= 1; --tmp2) {
            if (tmp1%tmp2 == 0) {
                break;
            }
        }
        nt[m] = tmp2;
        tmp1 /= nt[m];
        dx[m] = domain.getAxisSpan(m)/nt[m];
    }
#endif
    int i[domain.getNumDim()];
    memset(i, 0, sizeof(int)*domain.getNumDim()); i[0] = -1;
    for (t = tracers.begin(); t != tracers.end(); ++t) {
        // set spatial index
        i[0] += 1;
        for (int m = 0; m < domain.getNumDim()-1; ++m) {
            if (i[m] >= nt[m]) {
                i[m] = 0;
                i[m+1] += 1;
            }
        }
        // set coordinate
        for (int m = 0; m < domain.getNumDim(); ++m) {
            (*t)->getX(0)(m) = domain.getAxisStart(m)+dx[m]*0.5+dx[m]*i[m];
        }
        // set mesh index
        (*t)->getMeshIndex(0).locate(mesh, (*t)->getX(0));
    }
    // -------------------------------------------------------------------------
    REPORT_NOTICE(numTracer << " tracers are initialized.");
}

void TracerManager::output(const string &fileName, int timeLevel) {
    int ncId, numTracerDimId, numDimDimId;
    int qDimIds[2], qVarId;
    char str[100];
    int i;
    double *x;
    LADY_LIST<Tracer*>::iterator t;

    if (nc_create(fileName.c_str(), NC_CLOBBER, &ncId) != NC_NOERR) {
        REPORT_ERROR("Failed to open \"" << fileName << "\"!");
    }

    if (nc_def_dim(ncId, "num_tracer", tracers.size(), &numTracerDimId) != NC_NOERR) {
        REPORT_ERROR("Failed to define dimension \"num_tracer\"!");
    }

    if (nc_def_dim(ncId, "num_dim", domain->getNumDim(), &numDimDimId) != NC_NOERR) {

    }

    time_t curr_time;
    time(&curr_time);
    struct tm *timeinfo;
    timeinfo = gmtime(&curr_time);
    sprintf(str, "%4.2d-%2.2d-%2.2d",
            timeinfo->tm_year+1900,
            timeinfo->tm_mon+1,
            timeinfo->tm_mday+1);
    if (nc_put_att(ncId, NC_GLOBAL, "create_date", NC_CHAR, strlen(str), str)
        != NC_NOERR) {
        REPORT_ERROR("Failed to put attribute in \"" << fileName << "\"!");
    }

    qDimIds[0] = numTracerDimId; qDimIds[1] = numDimDimId;
    if (nc_def_var(ncId, "q", NC_DOUBLE, 2, qDimIds, &qVarId)
        != NC_NOERR) {
        REPORT_ERROR("Failed to define variable \"num_tracer\"!");
    }
    sprintf(str, "tracer centroid coordinates on %s", domain->getBrief().c_str());
    if (nc_put_att(ncId, qVarId, "long_name", NC_CHAR, strlen(str), str)
        != NC_NOERR) {
        REPORT_ERROR("Failed to put attribute in \"" << fileName << "\"!");
    }

    if (nc_enddef(ncId) != NC_NOERR) {

    }

    x = new double[tracers.size()*domain->getNumDim()];
    i = 0;
    for (t = tracers.begin(); t != tracers.end(); ++t) {
        for (int m = 0; m < domain->getNumDim(); ++m) {
            x[i++] = (*t)->getX(timeLevel)(m);
        }
    }
    if (nc_put_var(ncId, qVarId, x) != NC_NOERR) {
        REPORT_ERROR("Failed to put variable in \"" << fileName << "\"!");
    }

    if (nc_close(ncId) != NC_NOERR) {
        REPORT_ERROR("Failed to close file \"" << fileName << "\"!");
    }

    REPORT_NOTICE("File \"" << fileName << "\" is outputted.");
}

}
