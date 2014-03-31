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
protected:
    vector<double> density; //>! species density array
    TracerSkeleton *skeleton;

    /**
     *  Remapping parameters
     */
    int numConnectedCell;
    vector<TracerMeshCell*> connectedCells;
    double totalRemapWeight;
public:
    Tracer(int numDim);
    virtual ~Tracer();
    
    /**
     *  Add a species.
     */
    void addSpecies() { density.push_back(0); }

    double& getSpeciesDensity(int speciesIdx) {
#ifdef DEBUG
        if (speciesIdx >= density.size()) {
            REPORT_ERROR("Species index " << speciesIdx << " exceeds range [0," <<
                         density.size()-1 << "]!");
        }
#endif
        return density[speciesIdx];
    }

    double getSpeciesDensity(int speciesIdx) const {
#ifdef DEBUG
        if (speciesIdx >= density.size()) {
            REPORT_ERROR("Species index " << speciesIdx << " exceeds range [0," <<
                         density.size()-1 << "]!");
        }
#endif
        return density[speciesIdx];
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

    double getTotalRemapWeight() const {
        assert(totalRemapWeight != 0.0);
        return totalRemapWeight;
    }

    /**
     *  Update deformation matrix from tracer skeleton.
     *
     *  @param domain  the domain used for differencing coordinates.
     *  @param timeIdx the time level index.
     */
    void updateDeformMatrix(const LADY_DOMAIN &domain,
                            const LADY_MESH &mesh,
                            const TimeLevelIndex<2> &timeIdx);
    
    /**
     *  Check the linear deformation transformation validity.
     *
     *  @param domain  the domain.
     *  @param timeIdx the time level index.
     */
    void selfInspect(const LADY_DOMAIN &domain,
                     const TimeLevelIndex<2> &timeIdx);
};

}

#endif
