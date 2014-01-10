#include "TracerSkeleton.h"

namespace lady {

#define NUM_SKELETON_POINT numDim*2

TracerSkeleton::TracerSkeleton(Tracer *host, int numDim) {
    this->host = host;
    for (int l = 0; l < x.getNumLevel(); ++l) {
        x.getLevel(l).resize(NUM_SKELETON_POINT);
        idx.getLevel(l).resize(NUM_SKELETON_POINT);
        for (int i = 0; i < x.getLevel(l).size(); ++i) {
            x.getLevel(l)[i] = new LADY_SPACE_COORD(numDim);
            idx.getLevel(l)[i] = new LADY_MESH_INDEX(numDim);
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

void TracerSkeleton::init(const LADY_DOMAIN &domain, const LADY_MESH &mesh) {
    // set the body and initial spatial coordinates of skeleton points
    TimeLevelIndex<2> initTimeIdx;
    // -------------------------------------------------------------------------
    // sphere domain
    if (domain.getNumDim() == 2) {
        (*y[0])() << -0.5 <<  0.0 << arma::endr;
        (*y[1])() <<  0.5 <<  0.0 << arma::endr;
        (*y[2])() <<  0.0 << -0.5 << arma::endr;
        (*y[3])() <<  0.0 <<  0.5 << arma::endr;
        
        for (int i = 0; i < y.size(); ++i) {
            host->getSpaceCoord(domain, initTimeIdx, *y[i], *x.getLevel(0)[i]);
            idx.getLevel(0)[i]->locate(mesh, *x.getLevel(0)[i]);
        }
    } else if (domain.getNumDim() == 3) {
        REPORT_ERROR("Under construction!");
    }
}

vector<LADY_SPACE_COORD*>& TracerSkeleton::getXs(const TimeLevelIndex<2> &timeIdx) {
    return x.getLevel(timeIdx);
}

vector<LADY_BODY_COORD*>& TracerSkeleton::getYs() {
    return y;
}
    
vector<LADY_MESH_INDEX*>& TracerSkeleton::getIdxs(const TimeLevelIndex<2> &timeIdx) {
    return idx.getLevel(timeIdx);
}
    
}