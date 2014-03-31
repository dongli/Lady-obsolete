#ifndef __Lady_AdvectionManager__
#define __Lady_AdvectionManager__

#include "TracerManager.h"
#include "TracerMeshCell.h"

namespace lady {

// shortcuts for MLPACK classes
typedef mlpack::tree::BinarySpaceTree<
    mlpack::bound::HRectBound<2>,
    mlpack::range::RangeSearchStat> Tree;
typedef mlpack::metric::EuclideanDistance Metric;
typedef mlpack::range::RangeSearch<Metric, Tree> Searcher;

/**
 *  This class specifies the manager of linear advection.
 */
class AdvectionManager {
protected:
    TracerManager tracerManager;
    LADY_FIELD<TracerMeshCell> tracerMeshCells;
    LADY_REGRID *regrid; //>! used to interpolate velocity onto tracers
private:
    // range search parameters
    Tree *cellTree;                 //>! tree data structure for mesh cells for
                                    //>! avoiding rebuild of tree each time
    LADY_MATRIX cellCoords;         //>! collection of cell space coordinates
    vector<size_t> cellCoordsMap;   //>! mapping for cells since tree building
                                    //>! will modify the order of cells
public:
    AdvectionManager();
    ~AdvectionManager();

    /**
     *  Initialize advection manager.
     *
     *  @param domain     the spatial domain.
     *  @param mesh       the mesh.
     *  @param numParcelX the number of parcels (a.k.a., tracers) along x axis.
     *  @param numParcelY the number of parcels (a.k.a., tracers) along y axis.
     */
    void init(const LADY_DOMAIN &domain, const LADY_MESH &mesh,
              int numParcelX, int numParcelY);

    const LADY_FIELD<TracerMeshCell>& getTracerMeshCells() const {
        return tracerMeshCells;
    }

    /**
     *  Register a tracer species.
     *
     *  @param name  the name of the tracer.
     *  @param units the units of the tracer.
     *  @param brief the brief about the tracer.
     */
    void registerTracer(const string &name, const string &units,
                        const string &brief);

    /**
     *  Input one tracer species from given scalar field.
     *
     *  @param timeIdx the time level index.
     *  @param q       the input scalar field.
     */
    void input(const TimeLevelIndex<2> &timeIdx, vector<LADY_SCALAR_FIELD*> &q);

    /**
     *  Output tracers on old time level into netCDF file.
     *
     *  @param fileName   the file name.
     *  @param newTimeIdx the new time level index.
     */
    void output(const string &fileName, const TimeLevelIndex<2> &newTimeIdx);

    void diagnose(const TimeLevelIndex<2> &timeIdx);

    /**
     *  Advect tracers one time step forward, and remap tracers onto mesh cells.
     *
     *  @param dt         the time step size.
     *  @param oldTimeIdx the new time level index.
     *  @param V          the velocity field.
     */
    void advance(double dt, const TimeLevelIndex<2> &newTimeIdx,
                 const LADY_VELOCITY_FIELD &V);
private:
    /**
     *  Integrate the advection equations by using 4th-order Runge-Kutta method.
     *
     *  @param dt         the time step size.
     *  @param oldTimeIdx the old time level index.
     *  @param V          the velocity field.
     */
    void integrate_RK4(double dt, const TimeLevelIndex<2> &oldTimeIdx,
                       const LADY_VELOCITY_FIELD &V);

    /**
     *  Prepare the bidirectional remapping between tracers and mesh. Find out
     *  the mesh cells that a tracer will affected and calculate the weights.
     *
     *  @param timeIdx the time level index.
     */
    void connectTracersAndMesh(const TimeLevelIndex<2> &timeIdx);

    /**
     *  Remap the tracer mass from mesh cells to tracers.
     *
     *  @param timeIdx the time level index.
     */
    void remapMeshToTracers(const TimeLevelIndex<2> &timeIdx);
    
    /**
     *  Remap the tracer mass from tracers to mesh cells.
     *
     *  @param timeIdx the time level index.
     */
    void remapTracersToMesh(const TimeLevelIndex<2> &timeIdx);
};
}

#endif