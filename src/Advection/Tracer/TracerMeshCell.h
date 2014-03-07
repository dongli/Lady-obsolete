#ifndef __Lady_TracerMeshCell__
#define __Lady_TracerMeshCell__

#include "lady_commons.h"
#include "Tracer.h"

namespace lady {

class TracerMeshCell {
protected:
    LADY_SPACE_COORD *x; //>! center grid coordinate
    double volume; //>! cell volume for converting density and mass
    vector<double> m; //>! species mass array
    // remapping parameters
    int numConnectedTracer;
    vector<Tracer*> connectedTracers;
    vector<double> remapWeights;
    double totalRemapWeight;
public:
    TracerMeshCell();
    ~TracerMeshCell();

    /**
     *  Set the center coordinate.
     *
     *  @param x the grid coordinate.
     */
    void setCoord(const LADY_SPACE_COORD &x) { this->x = new LADY_SPACE_COORD(x); }

    /**
     *  Get the grid coordinate.
     *
     *  @return The grid coordinate.
     */
    const LADY_SPACE_COORD& getCoord() const { return *x; }

    /**
     *  Set the cell volume.
     *
     *  @param volume the cell volume.
     */
    void setVolume(double volume) { this->volume = volume; }

    /**
     *  Get the cell volume.
     *
     *  @return The cell volume.
     */
    double getVolume() const { return volume; }

    /**
     *  Add a species.
     */
    void addSpecies() { m.push_back(0.0); }

    double& getSpeciesMass(int s) { return m[s]; }
    double getSpeciesMass(int s) const { return m[s]; }

    /**
     *  Reset the connected tracer to empty for later updating.
     */
    void resetConnectedTracers();

    void resetSpeciesMass();

    void connect(Tracer *tracer, double weight);

    double getRemapWeight(Tracer *tracer) const;

    double getTotalRemapWeight() const { return totalRemapWeight; }
};

}

#endif