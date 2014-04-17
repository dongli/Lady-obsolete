#include "Tracer.h"
#include "TracerSkeleton.h"
#include "TracerMeshCell.h"

namespace lady {

Tracer::Tracer(int numDim) : Parcel(numDim) {
    skeleton = new TracerSkeleton(this, numDim);
    badType = GOOD_SHAPE;
    numConnectedCell = 0;
}

Tracer::~Tracer() {
    delete skeleton;
}

Tracer& Tracer::operator=(const Tracer &other) {
    Parcel::operator=(other);
    if (this != &other) {
        for (int s = 0; s < density.size(); ++s) {
            density[s] = other.density[s];
        }
        *skeleton = *(other.skeleton);
        // TODO: Handle connected cells.
    }
    return *this;
}

void Tracer::resetConnectedCells() {
    numConnectedCell = 0;
}

void Tracer::connect(TracerMeshCell *cell, double weight) {
    if (numConnectedCell == connectedCells.size()) {
        connectedCells.push_back(cell);
    } else {
        connectedCells[numConnectedCell] = cell;
    }
    numConnectedCell++;
}
    
void Tracer::updateDeformMatrix(const LADY_DOMAIN &domain,
                                const LADY_MESH &mesh,
                                const TimeLevelIndex<2> &timeIdx) {
    skeleton->updateLocalCoord(domain, timeIdx);
    const vector<vec> &xl = skeleton->getLocalCoords(timeIdx);
    const vector<LADY_BODY_COORD*> &y = skeleton->getBodyCoords();
    LADY_MATRIX &H0 = *H.getLevel(timeIdx);
    if (domain.getNumDim() == 2) {
        // elements of four matrices H_12, H_14, H_32, H_34
        double h11_1 = xl[0][0]/(*y[0])(0);
        double h21_1 = xl[0][1]/(*y[0])(0);
        double h11_3 = xl[2][0]/(*y[2])(0);
        double h21_3 = xl[2][1]/(*y[2])(0);
        double h12_2 = xl[1][0]/(*y[1])(1);
        double h22_2 = xl[1][1]/(*y[1])(1);
        double h12_4 = xl[3][0]/(*y[3])(1);
        double h22_4 = xl[3][1]/(*y[3])(1);
        // the final matrix is the combination of the four
        H0(0, 0) = (h11_1+h11_3)*0.5;
        H0(0, 1) = (h12_2+h12_4)*0.5;
        H0(1, 0) = (h21_1+h21_3)*0.5;
        H0(1, 1) = (h22_2+h22_4)*0.5;
        if (!svd(U, S, V, *H.getLevel(timeIdx))) {
            REPORT_ERROR("Encounter error with arma::svd!");
        }
#ifndef NDEBUG
        assert(S[0] >= S[1]);
#endif
        // reset the determinant of H to parcel volume
        // NOTE: Parcel volume is represented by the first tracer.
        if (density.size() != 0) {
            int s;
            for (s = 0; s < mass.size(); ++s) {
                if (mass[s] != 0) {
                    break;
                }
            }
            if (s != mass.size()) {
                double volume = mass[s]/density[s];
                double tmp = sqrt(volume/(S[0]*S[1]));
                S[0] *= tmp; S[1] *= tmp;
                H0 = U*diagmat(S)*V.t();
                resetSkeleton(domain, mesh, timeIdx);
            }
        }
        detH.getLevel(timeIdx) = S[0]*S[1];
        *invH.getLevel(timeIdx) = inv(H0);
        // check the degree of filamentation
        if (S[0]/S[1] > 100) {
            badType = EXTREME_FILAMENTATION;
        }
    } else if (domain.getNumDim() == 3) {
        REPORT_ERROR("Under construction!");
    }
    updateShapeSize(domain, timeIdx);
}

void Tracer::resetSkeleton(const LADY_DOMAIN &domain, const LADY_MESH &mesh,
                           const TimeLevelIndex<2> &timeIdx) {
    // reset skeleton points
    const vector<LADY_BODY_COORD*> &y = skeleton->getBodyCoords();
    vector<LADY_SPACE_COORD*> &x = skeleton->getSpaceCoords(timeIdx);
    vector<LADY_MESH_INDEX*> &meshIdx = skeleton->getMeshIdxs(timeIdx);
    for (int i = 0; i < x.size(); ++i) {
        getSpaceCoord(domain, timeIdx, *y[i], *x[i]);
        meshIdx[i]->reset();
        meshIdx[i]->locate(mesh, *x[i]);
    }
}

#ifndef NDEBUG
void Tracer::outputNeighbors(const TimeLevelIndex<2> &timeIdx,
                             const LADY_DOMAIN &domain) {
    std::ofstream file; file.open("neighbors.txt");
    // -------------------------------------------------------------------------
    // centroid location
    file << "centroid = (/" << getX(timeIdx)(0) << "," << getX(timeIdx)(1) << "/)" << endl;
    // -------------------------------------------------------------------------
    // neighbor cell locations
    if (numConnectedCell != 0) {
        file << "ngb_cells = new((/" << numConnectedCell << ",2/), double)" << endl;
        for (int m = 0; m < 2; ++m) {
            file << "ngb_cells(:," << m << ") = (/";
            for (int i = 0; i < numConnectedCell; ++i) {
                TracerMeshCell *cell = connectedCells[i];
                if (i != numConnectedCell-1) {
                    file << cell->getCoord()(m) << ",";
                } else {
                    file << cell->getCoord()(m) << "/)" << endl;
                }
            }
        }
    }
    // -------------------------------------------------------------------------
    // neighbor tracer locations
    int numNeighborTracer = 0;
    for (int i = 0; i < numConnectedCell; ++i) {
        numNeighborTracer += connectedCells[i]->getNumContainedTracer();
    }
    if (numNeighborTracer != 0) {
        file << "ngb_tracers = new((/" << numNeighborTracer << ",2/), double)" << endl;
        for (int m = 0; m < 2; ++m) {
            file << "ngb_tracers(:," << m << ") = (/";
            int k = 0;
            for (int i = 0; i < numConnectedCell; ++i) {
                TracerMeshCell *cell = connectedCells[i];
                vector<Tracer*> &tracers = cell->getContainedTracers();
                for (int j = 0; j < cell->getNumContainedTracer(); ++j) {
                    if (k != numNeighborTracer-1) {
                        file << tracers[j]->getX(timeIdx)(m) << ",";
                    } else {
                        file << tracers[j]->getX(timeIdx)(m) << "/)" << endl;
                    }
                    k++;
                }
            }
        }
    }
    // -------------------------------------------------------------------------
    // tracer shape
    int n = 100;
    file << "shape = new((/" << n << ",2/), double)" << endl;
    double dtheta = PI2/(n-1);
    LADY_BODY_COORD y(2);
    LADY_SPACE_COORD x(2);
    vector<vec::fixed<2> > shape(n);
    for (int i = 0; i < shape.size(); ++i) {
        double theta = i*dtheta;
        y(0) = cos(theta);
        y(1) = sin(theta);
        getSpaceCoord(domain, timeIdx, y, x);
        shape[i] = x();
    }
    for (int m = 0; m < 2; ++m) {
        file << "shape(:," << m << ") = (/";
        for (int i = 0; i < shape.size(); ++i) {
            if (i != shape.size()-1) {
                file << shape[i][m] << ",";
            } else {
                file << shape[i][m] << "/)" << endl;
            }
        }
    }
    file.close();
}
#endif

}
