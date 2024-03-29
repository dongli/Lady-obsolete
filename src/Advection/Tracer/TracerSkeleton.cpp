#include "TracerSkeleton.h"

namespace lady {

#define NUM_SKELETON_POINT numDim*2

TracerSkeleton::TracerSkeleton(Tracer *host, int numDim) {
    this->host = host;
    for (int l = 0; l < x.getNumLevel(); ++l) {
        x.getLevel(l).resize(NUM_SKELETON_POINT);
        idx.getLevel(l).resize(NUM_SKELETON_POINT);
        xl.getLevel(l).resize(NUM_SKELETON_POINT);
        for (int i = 0; i < x.getLevel(l).size(); ++i) {
            x.getLevel(l)[i] = new LADY_SPACE_COORD(numDim);
            idx.getLevel(l)[i] = new LADY_MESH_INDEX(numDim);
            xl.getLevel(l)[i].set_size(numDim);
        }
    }
    y.resize(NUM_SKELETON_POINT);
    for (int i = 0; i < y.size(); ++i) {
        y[i] = new LADY_BODY_COORD(numDim);
    }
}

TracerSkeleton::~TracerSkeleton() {
    for (int l = 0; l < x.getNumLevel(); ++l) {
        for (int i = 0; i < x.getLevel(l).size(); ++i) {
            delete x.getLevel(l)[i];
            delete idx.getLevel(l)[i];
        }
    }
    for (int i = 0; i < y.size(); ++i) {
        delete y[i];
    }
}
    
TracerSkeleton& TracerSkeleton::operator=(const TracerSkeleton &other) {
    if (this != &other) {
        for (int l = 0; l < x.getNumLevel(); ++l) {
            for (int i = 0; i < x.getLevel(l).size(); ++i) {
                *(x.getLevel(l)[i]) = *(other.x.getLevel(l)[i]);
                *(idx.getLevel(l)[i]) = *(other.idx.getLevel(l)[i]);
            }
        }
        for (int i = 0; i < y.size(); ++i) {
            *y[i] = *other.y[i];
        }
    }
    return *this;
}
    
/*
                  (0,1)
                    o
                    |
                    |
                    |
     (-1,0) o-------x------o (1,0)
                    |
                    |
                    |
                    o
                  (0,-1)
*/

void TracerSkeleton::init(const LADY_DOMAIN &domain, const LADY_MESH &mesh,
                          double size) {
    // set the body and initial spatial coordinates of skeleton points
    TimeLevelIndex<2> initTimeIdx;
    const LADY_SPACE_COORD &x0 = host->getX(initTimeIdx);
    if (domain.getNumDim() == 2) {
        double d = 0.5;
        (*y[0])() <<  -d << 0.0 << arma::endr;
        (*y[1])() << 0.0 <<  -d << arma::endr;
        (*y[2])() <<   d << 0.0 << arma::endr;
        (*y[3])() << 0.0 <<   d << arma::endr;
        double dtheta = PI2/y.size();
#ifdef LADY_USE_SPHERE_DOMAIN
        double lon, lat = M_PI_2-size/domain.getRadius();
        LADY_SPACE_COORD xr(domain.getNumDim());
        for (int i = 0; i < y.size(); ++i) {
            lon = i*dtheta;
            xr.setCoord(lon, lat);
            domain.rotateBack(x0, *x.getLevel(initTimeIdx)[i], xr);
            x.getLevel(initTimeIdx)[i]->transformToCart(domain); // TODO: Do we need this?
        }
#else
        REPORT_ERROR("Under construction!");
        for (int i = 0; i < y.size(); ++i) {
            double theta = i*dtheta;
            (*x.getLevel(initTimeIdx)[i])(0) = size*cos(theta)+x0(0);
            (*x.getLevel(initTimeIdx)[i])(1) = size*sin(theta)+x0(1);
        }
#endif
    } else if (domain.getNumDim() == 3) {
        REPORT_ERROR("Under construction!");
    }
    for (int i = 0; i < y.size(); ++i) {
        idx.getLevel(initTimeIdx)[i]->locate(mesh, *x.getLevel(initTimeIdx)[i]);
    }
}
    
void TracerSkeleton::updateLocalCoord(const LADY_DOMAIN &domain,
                                      const TimeLevelIndex<2> &timeIdx) {
    const LADY_SPACE_COORD &x0 = host->getX(timeIdx);
#ifdef LADY_USE_SPHERE_DOMAIN
    for (int i = 0; i < x.getLevel(timeIdx).size(); ++i) {
        domain.project(geomtk::SphereDomain::STEREOGRAPHIC, x0,
                       *x.getLevel(timeIdx)[i], xl.getLevel(timeIdx)[i]);
    }
#else
    for (int i = 0; i < x.getLevel(timeIdx).size(); ++i) {
        xl.getLevel(timeIdx)[i] = (*x.getLevel(timeIdx)[i])()-x0();
    }
#endif
}
    
}
