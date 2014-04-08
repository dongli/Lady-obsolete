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
    totalRemapWeight = 0.0;
}

void Tracer::connect(TracerMeshCell *cell, double weight) {
    if (numConnectedCell == connectedCells.size()) {
        connectedCells.push_back(cell);
    } else {
        connectedCells[numConnectedCell] = cell;
    }
    numConnectedCell++;
    totalRemapWeight += weight;
}
    
void Tracer::updateDeformMatrix(const LADY_DOMAIN &domain,
                                const LADY_MESH &mesh,
                                const TimeLevelIndex<2> &timeIdx) {
    skeleton->updateLocalCoord(domain, timeIdx);
    const vector<vec> &xl = skeleton->getLocalCoords(timeIdx);
    const vector<LADY_BODY_COORD*> &y = skeleton->getBodyCoords();
    const vector<LADY_SPACE_COORD*> &x = skeleton->getSpaceCoords(timeIdx);
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
        detH.getLevel(timeIdx) = det(H0);
        *invH.getLevel(timeIdx) = inv(H0);
        // check the degree of filamentation
        if (detH.getLevel(timeIdx) < 0) {
            // use last time step H
            if (!arma::svd(U, s, V, *H.getLevel(timeIdx-1))) {
                REPORT_ERROR("Encounter error with arma::svd!");
            }
            assert(s[0] >= s[1]);
            badType = EXTREME_FILAMENTATION;
            std::ofstream file("deform_matrix.dat");
            file << "H = new((/5,2,2/), double)" << endl;
            file << "H(0,:,:) = (/(/" << H0(0, 0) << "," << H0(0, 1) << "/), \\" << endl;
            file << "             (/" << H0(1, 0) << "," << H0(1, 1) << "/)/)" << endl;
            file << "H(1,:,:) = (/(/" << h11_1 << "," << h12_2 << "/), \\" << endl;
            file << "             (/" << h21_1 << "," << h22_2 << "/)/)" << endl;
            file << "H(2,:,:) = (/(/" << h11_1 << "," << h12_4 << "/), \\" << endl;
            file << "             (/" << h21_1 << "," << h22_4 << "/)/)" << endl;
            file << "H(3,:,:) = (/(/" << h11_3 << "," << h12_2 << "/), \\" << endl;
            file << "             (/" << h21_3 << "," << h22_2 << "/)/)" << endl;
            file << "H(4,:,:) = (/(/" << h11_3 << "," << h12_4 << "/), \\" << endl;
            file << "             (/" << h21_3 << "," << h22_4 << "/)/)" << endl;
            file << "xs = new((/4,2/), double)" << endl;
            file << "xs(0,:) = (/" << xl[0][0] << "," << xl[0][1] << "/)" << endl;
            file << "xs(1,:) = (/" << xl[1][0] << "," << xl[1][1] << "/)" << endl;
            file << "xs(2,:) = (/" << xl[2][0] << "," << xl[2][1] << "/)" << endl;
            file << "xs(3,:) = (/" << xl[3][0] << "," << xl[3][1] << "/)" << endl;
            file.close();
            CHECK_POINT;
        }
        // judge if the linear deformation approximation is ok
        if (badType == GOOD_SHAPE) {
            static const double threshold = 0.5;
            for (int i = 0; i < y.size(); ++i) {
                LADY_BODY_COORD Y(domain.getNumDim());
                getBodyCoord(domain, timeIdx, *x[i], Y);
                double bias = norm(Y()-(*y[i])(), 2);
                if (bias >= threshold) {
                    badType = POOR_LINEAR_APPROXIMATION;
                    std::ofstream file("deform_matrix.dat");
                    file << "H = new((/5,2,2/), double)" << endl;
                    file << "H(0,:,:) = (/(/" << H0(0, 0) << "," << H0(0, 1) << "/), \\" << endl;
                    file << "             (/" << H0(1, 0) << "," << H0(1, 1) << "/)/)" << endl;
                    file << "H(1,:,:) = (/(/" << h11_1 << "," << h12_2 << "/), \\" << endl;
                    file << "             (/" << h21_1 << "," << h22_2 << "/)/)" << endl;
                    file << "H(2,:,:) = (/(/" << h11_1 << "," << h12_4 << "/), \\" << endl;
                    file << "             (/" << h21_1 << "," << h22_4 << "/)/)" << endl;
                    file << "H(3,:,:) = (/(/" << h11_3 << "," << h12_2 << "/), \\" << endl;
                    file << "             (/" << h21_3 << "," << h22_2 << "/)/)" << endl;
                    file << "H(4,:,:) = (/(/" << h11_3 << "," << h12_4 << "/), \\" << endl;
                    file << "             (/" << h21_3 << "," << h22_4 << "/)/)" << endl;
                    file << "xs = new((/4,2/), double)" << endl;
                    file << "xs(0,:) = (/" << xl[0][0] << "," << xl[0][1] << "/)" << endl;
                    file << "xs(1,:) = (/" << xl[1][0] << "," << xl[1][1] << "/)" << endl;
                    file << "xs(2,:) = (/" << xl[2][0] << "," << xl[2][1] << "/)" << endl;
                    file << "xs(3,:) = (/" << xl[3][0] << "," << xl[3][1] << "/)" << endl;
                    file.close();
                    CHECK_POINT;
                    break;
                }
            }
        }
    } else if (domain.getNumDim() == 3) {
        REPORT_ERROR("Under construction!");
    }
    updateShapeSize(domain, timeIdx);
}

void Tracer::mixWithNeighborTracers(const TimeLevelIndex<2> &timeIdx,
                                    const LADY_DOMAIN &domain) {
    vector<Tracer*> tracers;
    for (int i = 0; i < numConnectedCell; ++i) {
        TracerMeshCell *cell = connectedCells[i];
        for (int j = 0; j < cell->getNumContainedTracer(); ++j) {
            Tracer *ngb = cell->getContainedTracers()[j];
            tracers.push_back(ngb);
        }
    }
    LADY_BODY_COORD y(domain.getNumDim());
    vec weights(tracers.size());
    mat newDensity(density.size(), tracers.size());
    for (int i = 0; i < tracers.size(); ++i) {
        for (int s = 0; s < density.size(); ++s) {
            newDensity(s, i) = 0;
        }
        for (int j = 0; j < tracers.size(); ++j) {
            tracers[i]->getBodyCoord(domain, timeIdx, tracers[j]->getX(timeIdx), y);
            weights[j] = tracers[i]->getShapeFunction(timeIdx, y);
        }
        weights = weights/sum(weights);
        for (int j = 0; j < tracers.size(); ++j) {
            for (int s = 0; s < density.size(); ++s) {
                newDensity(s, i) += tracers[j]->density[s]*weights[j];
            }
        }
    }
    for (int i = 0; i < tracers.size(); ++i) {
        for (int s = 0; s < density.size(); ++s) {
            tracers[i]->density[s] = newDensity(s, i);
        }
    }
    CHECK_POINT;
}

void Tracer::adjustDeformMatrix(const TimeLevelIndex<2> &timeIdx) {
    assert(badType == EXTREME_FILAMENTATION);
    // reduce the long and short semiaxis ratio by half times
    s[0] *= 0.6;
    s[1] *= 1.2;
    (*H.getLevel(timeIdx)) = U*diagmat(s)*V.t();
    detH.getLevel(timeIdx) = det(*H.getLevel(timeIdx));
    *invH.getLevel(timeIdx) = inv(*H.getLevel(timeIdx));
    assert(detH.getLevel(timeIdx) > 0);
}

void Tracer::resetSkeleton(const LADY_DOMAIN &domain, const LADY_MESH &mesh,
                           const TimeLevelIndex<2> &timeIdx) {
    assert(badType != GOOD_SHAPE);
    // reset skeleton points
    const vector<LADY_BODY_COORD*> &y = skeleton->getBodyCoords();
    vector<LADY_SPACE_COORD*> &x = skeleton->getSpaceCoords(timeIdx);
    vector<LADY_MESH_INDEX*> &meshIdx = skeleton->getMeshIdxs(timeIdx);
    for (int i = 0; i < x.size(); ++i) {
        getSpaceCoord(domain, timeIdx, *y[i], *x[i]);
        meshIdx[i]->locate(mesh, *x[i]);
    }
    badType = GOOD_SHAPE;
}

#ifdef DEBUG
void Tracer::outputNeighbors(const TimeLevelIndex<2> &timeIdx,
                             const LADY_DOMAIN &domain) {
    std::ofstream file; file.open("neighbors.txt");
    // -------------------------------------------------------------------------
    // centroid location
    file << "centroid = (/" << getX(timeIdx)(0) << "," << getX(timeIdx)(1) << "/)" << endl;
    // -------------------------------------------------------------------------
    // neighbor cell locations
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
    // -------------------------------------------------------------------------
    // neighbor tracer locations
    int numNeighborTracer = 0;
    for (int i = 0; i < numConnectedCell; ++i) {
        numNeighborTracer += connectedCells[i]->getNumContainedTracer();
    }
    
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
        y(0) = 1.2*cos(theta);
        y(1) = 1.2*sin(theta);
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
