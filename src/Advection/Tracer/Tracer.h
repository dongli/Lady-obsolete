#ifndef __Lady_Tracer__
#define __Lady_Tracer__

#include "lady_commons.h"

namespace lady {

class TracerSkeleton;

/**
 *  This class describes the tracer that is used to be advected by external wind
 *  field.
 */
class Tracer {
protected:
    int ID;
    TimeLevels<LADY_SPACE_COORD*, 2> q; // centroid coordinate
    TimeLevels<LADY_MATRIX*, 2> H; // deformation matrix
    TimeLevels<LADY_MESH_INDEX*, 2> idx; // tracer mesh index
    TimeLevels<vec, 2> ms; // species mass array

    TimeLevels<double, 2> detH; // determinant of H
    TimeLevels<LADY_MATRIX*, 2> invH; // inversion of H

    TracerSkeleton *skeleton;
public:
    Tracer(int numDim);
    virtual ~Tracer();
    
    /**
     *  Add a species.
     */
    void addSpecies();

    Tracer& operator=(const Tracer &other);
    
    /**
     *  Set the ID of tracer.
     *
     *  @param ID the given ID.
     */
    void setID(int ID);

    /**
     *  Get the ID of tracer.
     *
     *  @return The ID.
     */
    int getID() const;

    /**
     *  Get the spatial coordinate of tracer.
     *
     *  @param timeIdx the time level index.
     *
     *  @return The spatial coordinate.
     */
    LADY_SPACE_COORD& getX(const TimeLevelIndex<2> &timeIdx);

    /**
     *  Get the deformation matrix of tracer.
     *
     *  @param timeIdx the time level index.
     *
     *  @return The deformation matrix.
     */
    LADY_MATRIX& getH(const TimeLevelIndex<2> &timeIdx);

    /**
     *  Get the inverse deformation matrix of tracer.
     *
     *  @param timeIdx the time level index.
     *
     *  @return The inverse deformation matrix.
     */
    LADY_MATRIX& getInvH(const TimeLevelIndex<2> &timeIdx);

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

    /**
     *  Transform body coordinate into spatial coordinate.
     *
     *  @param domain  the domain used for differencing coordinates.
     *  @param timeIdx the time level index.
     *  @param y       the body coordinate.
     *  @param x       the spatial coordinate.
     */
    void getSpaceCoord(const LADY_DOMAIN &domain,
                       const TimeLevelIndex<2> &timeIdx,
                       const LADY_BODY_COORD &y, LADY_SPACE_COORD &x);
    
    /**
     *  Transform spatial coordinate into body coordinate.
     *
     *  @param domain  the domain used for differencing coordinates.
     *  @param timeIdx the time level index.
     *  @param x       the spatial coordinate.
     *  @param y       the body coordinate.
     */
    void getBodyCoord(const LADY_DOMAIN &domain,
                      const TimeLevelIndex<2> &timeIdx,
                      const LADY_SPACE_COORD &x, LADY_BODY_COORD &y);

    double getShapeFunction(const TimeLevelIndex<2> &timeIdx,
                            const LADY_BODY_COORD &y);

    /**
     *  Update deformation matrix from tracer skeleton.
     *
     *  @param domain  the domain used for differencing coordinates.
     *  @param timeIdx the time level index.
     */
    void updateH(const LADY_DOMAIN &domain,
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
