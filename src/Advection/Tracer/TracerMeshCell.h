#ifndef __Lady_TracerMeshCell__
#define __Lady_TracerMeshCell__

#include "lady_commons.h"
#include "Tracer.h"

namespace lady {

class TracerMeshCell {
protected:
    LADY_SPACE_COORD *x; //>! center grid coordinate
    double volume; //>! cell volume for converting density and mass
    vector<double> ms; //>! species mass array
    unordered_map<Tracer*, double> tracers; //>! tracer-weight map
    double totalRemapWeight;
public:
    TracerMeshCell();
    ~TracerMeshCell();

    /**
     *  Set the center coordinate.
     *
     *  @param x the grid coordinate.
     */
    void setCoord(const LADY_SPACE_COORD &x);

    /**
     *  Get the grid coordinate.
     *
     *  @return The grid coordinate.
     */
    const LADY_SPACE_COORD& getCoord() const;

    /**
     *  Set the cell volume.
     *
     *  @param volume the cell volume.
     */
    void setVolume(double volume);

    /**
     *  Get the cell volume.
     *
     *  @return The cell volume.
     */
    double getVolume() const;

    /**
     *  Add a species.
     */
    void addSpecies();

    double& getSpeciesMass(int speciesIdx);
    double getSpeciesMass(int speciesIdx) const;

    void resetConnectedTracers();

    void connect(Tracer *tracer, double weight);

    double getWeight(Tracer *tracer) const;

    double getTotalRemapWeight() const;
};

}

#endif