#ifndef __Lady_Tracer__
#define __Lady_Tracer__

#include "lady_commons.h"
#include "Parcel.h"

namespace lady {

class TracerSkeleton;
class TracerMeshCell;

/**
 *  This class describes the tracer that is used to be advected by external wind
 *  field. It derives from parcel.
 */
class Tracer : public Parcel {
public:
    enum BadDeformMatrixType {
        GOOD_SHAPE, EXTREME_FILAMENTATION, NOT_RESOLVED
    };
protected:
    vector<double> density; //>! species density array
    vector<double> mass;    //>! species mass array
    TracerSkeleton *skeleton;
    BadDeformMatrixType badType;
    /**
     *  Remapping parameters
     */
    int numConnectedCell;
    vector<TracerMeshCell*> connectedCells;

    TracerMeshCell *hostCell;
public:
    Tracer(int numDim);
    virtual ~Tracer();
    
    /**
     *  Add a species.
     */
    void addSpecies() { density.push_back(0); mass.push_back(0); }

    double& getSpeciesDensity(int speciesIdx) {
#ifndef NDEBUG
        if (speciesIdx >= density.size()) {
            REPORT_ERROR("Species index " << speciesIdx << " exceeds range [0," <<
                         density.size()-1 << "]!");
        }
#endif
        return density[speciesIdx];
    }

    double getSpeciesDensity(int speciesIdx) const {
#ifndef NDEBUG
        if (speciesIdx >= density.size()) {
            REPORT_ERROR("Species index " << speciesIdx << " exceeds range [0," <<
                         density.size()-1 << "]!");
        }
#endif
        return density[speciesIdx];
    }
    
    void calcSpeciesMass(const TimeLevelIndex<2> &timeIdx, int speciesIdx) {
        mass[speciesIdx] = density[speciesIdx]*detH.getLevel(timeIdx);
    }

    double& getSpeciesMass(int speciesIdx) {
        return mass[speciesIdx];
    }

    double getSpeciesMass(int speciesIdx) const {
        return mass[speciesIdx];
    }

    void resetSpecies() {
        for (int s = 0; s < density.size(); ++s) {
            density[s] = 0.0;
        }
    }

    Tracer& operator=(const Tracer &other);

    /**
     *  Get the mesh index of tracer.
     *
     *  @param timeIdx the time level index.
     *
     *  @return The mesh index.
     */
    LADY_MESH_INDEX& getMeshIndex(const TimeLevelIndex<2> &timeIdx) {
        return *idx.getLevel(timeIdx);
    }

    /**
     *  Get the skeleton of tracer.
     *
     *  @return The skeleton object.
     */
    TracerSkeleton& getSkeleton() { return *skeleton; }

    void setBadType(BadDeformMatrixType badType) { this->badType = badType; }

    BadDeformMatrixType getBadType() const { return badType; }

    /**
     *  Reset the connected cells to empty for later updating.
     */
    void resetConnectedCells();

    /**
     *  Connect the tracer with the given mesh cell.
     *
     *  @param cell   the mesh cell to be connected.
     *  @param weight the remapping weight.
     */
    void connect(TracerMeshCell *cell, double weight);

    /**
     *  Get the connected mesh cells.
     *
     *  @return The connected mesh cell list.
     */
    vector<TracerMeshCell*>& getConnectedCells() { return connectedCells; }
    
    int getNumConnectedCell() const { return numConnectedCell; }

    void setHostCell(TracerMeshCell *cell) { hostCell = cell; }
    
    TracerMeshCell* getHostCell() const { return hostCell; }

    /**
     *  Update deformation matrix from tracer skeleton.
     *
     *  @param domain  the spatial domain.
     *  @param mesh    the model mesh.
     *  @param timeIdx the time level index.
     */
    void updateDeformMatrix(const LADY_DOMAIN &domain,
                            const LADY_MESH &mesh,
                            const TimeLevelIndex<2> &timeIdx);

    /**
     *  Reset deformation matrix from given major axis vertex and scale matrix.
     *
     *  @param domain  the spatial domain.
     *  @param mesh    the model mesh.
     *  @param timeIdx the time level index.
     *  @param x       the major axis vertex spatial coordinate.
     *  @param S       the scale matrix.
     */
    void resetDeformMatrix(const LADY_DOMAIN &domain,
                           const LADY_MESH &mesh,
                           const TimeLevelIndex<2> &timeIdx,
                           const LADY_SPACE_COORD &x, const vec &S);

    /**
     *  Reset tracer skeleton points from deformation matrix
     *
     *  @param domain  the spatial domain.
     *  @param mesh    the model mesh.
     *  @param timeIdx the time level index.
     */
    void resetSkeleton(const LADY_DOMAIN &domain, const LADY_MESH &mesh,
                       const TimeLevelIndex<2> &timeIdx);

#ifndef NDEBUG
    void outputNeighbors(const TimeLevelIndex<2> &timeIdx,
                         const LADY_DOMAIN &domain);
#endif
};

}

#endif
