#ifndef __Lady_TracerMeshCell__
#define __Lady_TracerMeshCell__

#include "lady_commons.h"
#include "Tracer.h"

namespace lady {

class TracerMeshCell {
protected:
    int ID;
    LADY_SPACE_COORD *x;    //>! center grid coordinate
    double volume;          //>! cell volume for converting density and mass
    vector<double> density; //>! species density array
    // remapping parameters
    int numConnectedTracer;
    vector<Tracer*> connectedTracers;
    vector<double> remapWeights;
    // mixing parameters
    int numContainedTracer;
    vector<Tracer*> containedTracers;
public:
    TracerMeshCell();
    ~TracerMeshCell();

    void setID(int ID) { this->ID = ID; }
    
    int getID() const { return ID; }

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
    void addSpecies() { density.push_back(0); }

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

    void resetSpecies() {
        for (int s = 0; s < density.size(); ++s) {
            density[s] = 0;
        }
    }

    /**
     *  Reset the connected tracer to empty for later updating.
     */
    void resetConnectedTracers();

    void connect(Tracer *tracer, double weight);

    void disconnect(Tracer *tracer);

    int getNumConnectedTracer() const { return numConnectedTracer; }

    vector<Tracer*>& getConnectedTracers() { return connectedTracers; }

    double getRemapWeight(Tracer *tracer) const;

    void resetContainedTracers();

    void contain(Tracer *tracer);

    void discontain(Tracer *tracer);

    int getNumContainedTracer() const { return numContainedTracer; }

    vector<Tracer*>& getContainedTracers() { return containedTracers; }
};

}

#endif