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
    TimeLevels<vector<double>, 2> ms; //>! species mass array
    TracerSkeleton *skeleton;

    /**
     *  Remapping parameters
     */
    list<TracerMeshCell*> cells;
    double totalRemapWeight;
public:
    Tracer(int numDim);
    virtual ~Tracer();
    
    /**
     *  Add a species.
     */
    void addSpecies();

    double& getSpeciesMass(const TimeLevelIndex<2> &timeIdx, int speciesIdx);
    double getSpeciesMass(const TimeLevelIndex<2> &timeIdx, int speciesIdx) const;

    Tracer& operator=(const Tracer &other);

    /**
     *  Get the mesh index of tracer.
     *
     *  @param timeIdx the time level index.
     *
     *  @return The mesh index.
     */
    LADY_MESH_INDEX& getMeshIndex(const TimeLevelIndex<2> &timeIdx);

    /**
     *  Get the skeleton of tracer.
     *
     *  @return The skeleton object.
     */
    TracerSkeleton& getSkeleton();

    void resetConnectedCells();

    /**
     *  Connect the tracer with the given mesh cell.
     *
     *  @param cell the mesh cell to be connected.
     */
    void connect(TracerMeshCell *cell);

    /**
     *  Get the connected mesh cells.
     *
     *  @return The connected mesh cell list.
     */
    const list<TracerMeshCell*>& getConnectedCells() const;

    double getTotalRemapWeight() const;

    /**
     *  Update deformation matrix from tracer skeleton.
     *
     *  @param domain  the domain used for differencing coordinates.
     *  @param timeIdx the time level index.
     */
    void updateDeformMatrix(const LADY_DOMAIN &domain,
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
