#include <TracerManager.h>

namespace lady
{

TracerManager::TracerManager() {
    REPORT_ONLINE;
}

TracerManager::~TracerManager() {
    REPORT_OFFLINE;
}

void TracerManager::init(const LADY_DOMAIN &domain, const LADY_MESH &mesh,
                         int numTracer) {
    this->domain = &domain;
    tracers.resize(numTracer);
    int ID = 0;
    LADY_LIST<Tracer*>::iterator tracer = tracers.begin();
    for (; tracer != tracers.end(); ++tracer) {
        *tracer = new Tracer(domain.getNumDim());
        (*tracer)->setID(ID++);
    }
    // -------------------------------------------------------------------------
    // set initial centroid coordinates and deformation matrices of tracers
    vec dx(domain.getNumDim());
    int nt[domain.getNumDim()];
#define DEBUG_CROSS_POLE 0
#define DEBUG_EQUATOR 1
#if DEBUG_CROSS_POLE == 1
    assert(dynamic_cast<const geomtk::SphereDomain*>(&domain));
    nt[0] = numTracer; nt[1] = 1;
    dx[0] = PI2/numTracer; dx[1] = 3*RAD;
#elif DEBUG_EQUATOR == 1
    assert(dynamic_cast<const geomtk::SphereDomain*>(&domain));
    nt[0] = numTracer; nt[1] = 1;
    dx[0] = PI2/numTracer; dx[1] = 180*RAD;
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
    for (tracer = tracers.begin(); tracer != tracers.end(); ++tracer) {
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
            (*tracer)->getX(0)(m) = domain.getAxisStart(m)+dx[m]*0.5+dx[m]*i[m];
        }
        // set mesh index
        (*tracer)->getMeshIndex(0).locate(mesh, (*tracer)->getX(0));
        // set deformation matrix
        (*tracer)->getH(0).zeros();
#if DEBUG_CROSS_POLE == 1 || DEBUG_EQUATOR == 1
        (*tracer)->getH(0)(0, 0) = 2*dx[0];
        (*tracer)->getH(0)(1, 1) = 2*dx[0];
#else
        for (int m = 0; m < domain.getNumDim(); ++m) {
            (*tracer)->getH(0)(m, m) = 2*dx.max(); // TODO: Could we set H this way?
        }
#endif
    }
    // -------------------------------------------------------------------------
    REPORT_NOTICE(numTracer << " tracers are initialized.");
}

void TracerManager::output(const string &fileName, int timeLevel) {
    int ncId, numTracerDimId, numDimDimId;
    int qDimIds[2], qVarId;
    int hDimIds[3], hVarId;
    char str[100];
    int i;
    double *x;
    LADY_LIST<Tracer*>::iterator tracer;

    if (nc_create(fileName.c_str(), NC_CLOBBER, &ncId) != NC_NOERR) {
        REPORT_ERROR("Failed to open \"" << fileName << "\"!");
    }

    if (nc_def_dim(ncId, "num_tracer", tracers.size(), &numTracerDimId)
        != NC_NOERR) {
        REPORT_ERROR("Failed to define dimension \"num_tracer\"!");
    }

    if (nc_def_dim(ncId, "num_dim", domain->getNumDim(), &numDimDimId)
        != NC_NOERR) {

    }

    time_t curr_time;
    time(&curr_time);
    struct tm *timeinfo;
    timeinfo = gmtime(&curr_time);
    sprintf(str, "UTC %4.2d-%2.2d-%2.2d",
            timeinfo->tm_year+1900,
            timeinfo->tm_mon+1,
            timeinfo->tm_mday);
    if (nc_put_att(ncId, NC_GLOBAL, "create_date", NC_CHAR, strlen(str), str)
        != NC_NOERR) {
        REPORT_ERROR("Failed to put attribute in \"" << fileName << "\"!");
    }

    qDimIds[0] = numTracerDimId;
    qDimIds[1] = numDimDimId;
    if (nc_def_var(ncId, "q", NC_DOUBLE, 2, qDimIds, &qVarId) != NC_NOERR) {
        REPORT_ERROR("Failed to define variable \"q\"!");
    }
    sprintf(str, "tracer centroid coordinates on %s", domain->getBrief().c_str());
    if (nc_put_att(ncId, qVarId, "long_name", NC_CHAR, strlen(str), str)
        != NC_NOERR) {
        REPORT_ERROR("Failed to put attribute in \"" << fileName << "\"!");
    }

    hDimIds[0] = numTracerDimId;
    hDimIds[1] = numDimDimId;
    hDimIds[2] = numDimDimId;
    if (nc_def_var(ncId, "h", NC_DOUBLE, 3, hDimIds, &hVarId) != NC_NOERR) {
        REPORT_ERROR("Failed to define variable \"q\"!");
    }
    sprintf(str, "tracer deformation matrix");
    if (nc_put_att(ncId, hVarId, "long_name", NC_CHAR, strlen(str), str)
        != NC_NOERR) {
        REPORT_ERROR("Failed to put attribute in \"" << fileName << "\"!");
    }
    
    if (nc_enddef(ncId) != NC_NOERR) {

    }
    // -------------------------------------------------------------------------
    x = new double[tracers.size()*domain->getNumDim()];
    i = 0;
    for (tracer = tracers.begin(); tracer != tracers.end(); ++tracer) {
        for (int m = 0; m < domain->getNumDim(); ++m) {
            x[i++] = (*tracer)->getX(timeLevel)(m);
        }
    }
    if (nc_put_var(ncId, qVarId, x) != NC_NOERR) {
        REPORT_ERROR("Failed to put variable in \"" << fileName << "\"!");
    }

    x = new double[tracers.size()*domain->getNumDim()*domain->getNumDim()];
    i = 0;
    for (tracer = tracers.begin(); tracer != tracers.end(); ++tracer) {
        for (int m1 = 0; m1 < domain->getNumDim(); ++m1) {
            for (int m2 = 0; m2 < domain->getNumDim(); ++m2) {
                x[i++] = (*tracer)->getH(timeLevel)(m1, m2);
            }
        }
    }
    if (nc_put_var(ncId, hVarId, x) != NC_NOERR) {
        REPORT_ERROR("Failed to put variable in \"" << fileName << "\"!");
    }
    // -------------------------------------------------------------------------
    if (nc_close(ncId) != NC_NOERR) {
        REPORT_ERROR("Failed to close file \"" << fileName << "\"!");
    }

    REPORT_NOTICE("File \"" << fileName << "\" is outputted.");
}

}
