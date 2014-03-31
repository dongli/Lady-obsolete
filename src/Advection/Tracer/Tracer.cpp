#include "Tracer.h"
#include "TracerSkeleton.h"
#include "TracerMeshCell.h"

namespace lady {

Tracer::Tracer(int numDim) : Parcel(numDim) {
    for (int l = 0; l < idx.getNumLevel(); ++l) {
        idx.getLevel(l) = new LADY_MESH_INDEX(numDim);
    }
    skeleton = new TracerSkeleton(this, numDim);
    numConnectedCell = 0;
}

Tracer::~Tracer() {
    for (int l = 0; l < idx.getNumLevel(); ++l) {
        delete idx.getLevel(l);
    }
    delete skeleton;
}
    
void Tracer::addSpecies() {
    m.push_back(0.0);
}

double& Tracer::getSpeciesMass(int s) {
#ifdef DEBUG
    if (s >= m.size()) {
        REPORT_ERROR("Species index " << s << " exceeds range [0," <<
                     m.size()-1 << "]!");
    }
#endif
    return m[s];
}

double Tracer::getSpeciesMass(int s) const {
#ifdef DEBUG
    if (s >= m.size()) {
        REPORT_ERROR("Species index " << s << " exceeds range [0," <<
                     m.size()-1 << "]!");
    }
#endif
    return m[s];
}

Tracer& Tracer::operator=(const Tracer &other) {
    Parcel::operator=(other);
    if (this != &other) {
        for (int l = 0; l < idx.getNumLevel(); ++l) {
            *(idx.getLevel(l)) = *(other.idx.getLevel(l));
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

void Tracer::resetSpeciesMass() {
    for (int s = 0; s < m.size(); ++s) {
        m[s] = 0.0;
    }
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

double Tracer::getTotalRemapWeight() const {
    assert(totalRemapWeight != 0.0);
    return totalRemapWeight;
}
    
void Tracer::updateDeformMatrix(const LADY_DOMAIN &domain,
                                const LADY_MESH &mesh,
                                const TimeLevelIndex<2> &timeIdx) {
    skeleton->updateLocalCoord(domain, timeIdx);
    const vector<vec> &xl = skeleton->getLocalCoords(timeIdx);
    const vector<LADY_BODY_COORD*> &y = skeleton->getBodyCoords();
    LADY_MATRIX &H0 = *H.getLevel(timeIdx);
    if (domain.getNumDim() == 2) {
        // weights for the four matrix members
        double w1 = 0.5;
        double w2 = 0.5;
        double w3 = 0.5;
        double w4 = 0.5;
        // elements of four matrices H_12, H_14, H_32, H_34
        double h11_1 = xl[0][0]/(*y[0])(0);
        double h21_1 = xl[0][1]/(*y[0])(0);
        double h11_3 = xl[2][0]/(*y[2])(0);
        double h21_3 = xl[2][1]/(*y[2])(0);
        double h12_2 = xl[1][0]/(*y[1])(1);
        double h22_2 = xl[1][1]/(*y[1])(1);
        double h12_4 = xl[3][0]/(*y[3])(1);
        double h22_4 = xl[3][1]/(*y[3])(1);
        
        while (true) {
            H0(0, 0) = (h11_1*w1+h11_3*w3)/(w1+w3);
            H0(0, 1) = (h12_2*w2+h12_4*w4)/(w2+w4);
            H0(1, 0) = (h21_1*w1+h21_3*w3)/(w1+w3);
            H0(1, 1) = (h22_2*w2+h22_4*w4)/(w2+w4);
            detH.getLevel(timeIdx) = det(H0);
            if (detH.getLevel(timeIdx) < 0) {
                double detH_12 = h11_1*h22_2-h21_1*h12_2;
                double detH_14 = h11_1*h22_4-h21_1*h12_4;
                double detH_32 = h11_3*h22_2-h21_3*h12_2;
                double detH_34 = h11_3*h22_4-h21_3*h12_4;
                if (detH_12 < 0 && detH_14 < 0 && detH_32 < 0 && detH_34 < 0) {
                    REPORT_ERROR("Tracer " << ID << " has bad deformation " <<
                                 "matrix with negative determinant, and can " <<
                                 "not be adjusted!");
                }
                if (detH_12 < 0 || detH_14 < 0) {
                    w1 *= 0.5;
                }
                if (detH_12 < 0 || detH_32 < 0) {
                    w2 *= 0.5;
                }
                if (detH_32 < 0 || detH_34 < 0) {
                    w3 *= 0.5;
                }
                if (detH_14 < 0 || detH_34 < 0) {
                    w4 *= 0.5;
                }
                REPORT_WARNING("Tracer " << ID << " encounter bad deformation " <<
                               "matrix, but can be adjusted.");
            } else {
                break;
            }
        }
        
//        if (ID == 1963) {
//            std::ofstream file("deform_matrix.dat");
//            file << "H = new((/5,2,2/), double)" << endl;
//            file << "H(0,:,:) = (/(/" << H0(0, 0) << "," << H0(0, 1) << "/), \\" << endl;
//            file << "             (/" << H0(1, 0) << "," << H0(1, 1) << "/)/)" << endl;
//            file << "H(1,:,:) = (/(/" << h11_1 << "," << h12_2 << "/), \\" << endl;
//            file << "             (/" << h21_1 << "," << h22_2 << "/)/)" << endl;
//            file << "H(2,:,:) = (/(/" << h11_1 << "," << h12_4 << "/), \\" << endl;
//            file << "             (/" << h21_1 << "," << h22_4 << "/)/)" << endl;
//            file << "H(3,:,:) = (/(/" << h11_3 << "," << h12_2 << "/), \\" << endl;
//            file << "             (/" << h21_3 << "," << h22_2 << "/)/)" << endl;
//            file << "H(4,:,:) = (/(/" << h11_3 << "," << h12_4 << "/), \\" << endl;
//            file << "             (/" << h21_3 << "," << h22_4 << "/)/)" << endl;
//            file << "xs = new((/4,2/), double)" << endl;
//            file << "xs(0,:) = (/" << xl[0][0] << "," << xl[0][1] << "/)" << endl;
//            file << "xs(1,:) = (/" << xl[1][0] << "," << xl[1][1] << "/)" << endl;
//            file << "xs(2,:) = (/" << xl[2][0] << "," << xl[2][1] << "/)" << endl;
//            file << "xs(3,:) = (/" << xl[3][0] << "," << xl[3][1] << "/)" << endl;
//            file.close();
//            CHECK_POINT
//        }
    } else if (domain.getNumDim() == 3) {
        REPORT_ERROR("Under construction!");
    }
    *invH.getLevel(timeIdx) = inv(H0);
    updateShapeSize(domain, timeIdx);

//    if (ID == 115) {
//        std::ofstream file("y.dat");
//        for (int i = 0; i < x.size(); ++i) {
//            LADY_BODY_COORD y(2);
//            getBodyCoord(domain, timeIdx, *skeleton->getSpaceCoords(timeIdx)[i], y);
//            file << y() << endl;
//        }
//        file.close();
//        CHECK_POINT
//    }
    
//    // reset skeleton points
//    vector<LADY_SPACE_COORD*> &x = skeleton->getSpaceCoords(timeIdx);
//    const vector<LADY_BODY_COORD*> &y = skeleton->getBodyCoords();
//    vector<LADY_MESH_INDEX*> &meshIdx = skeleton->getMeshIdxs(timeIdx);
//    for (int i = 0; i < x.size(); ++i) {
//        getSpaceCoord(domain, timeIdx, *y[i], *x[i]);
//        meshIdx[i]->locate(mesh, *x[i]);
//    }
}

void Tracer::selfInspect(const LADY_DOMAIN &domain,
                         const TimeLevelIndex<2> &timeIdx) {
    static const double threshold = 0.2;
    const vector<LADY_BODY_COORD*> &ys = skeleton->getBodyCoords();
    vector<LADY_SPACE_COORD*> &xs = skeleton->getSpaceCoords(timeIdx);
    for (int i = 0; i < ys.size(); ++i) {
        LADY_BODY_COORD y(domain.getNumDim());
        getBodyCoord(domain, timeIdx, *xs[i], y);
        double bias = norm(y()-(*ys[i])(), 2);
        if (bias >= threshold) {
            cout << bias << endl;
            ys[i]->print();
            y.print();
            cout << ID << endl;
        }
    }
}

}
