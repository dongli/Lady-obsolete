#ifndef __Lady_Parcel__
#define __Lady_Parcel__

#include "lady_commons.h"

namespace lady {

/**
 *  This class specifies the parcel, which is the computational unit of Lady.
 */
class Parcel {
protected:
    int ID;
    TimeLevels<LADY_SPACE_COORD*, 2> q;     //>! centroid coordinate
    TimeLevels<LADY_MATRIX*, 2> H;          //>! deformation matrix
    TimeLevels<double, 2> detH;             //>! determinant of H
    TimeLevels<LADY_MATRIX*, 2> invH;       //>! inversion of H
    TimeLevels<vec::fixed<2>, 2> shapeSize; //>! parcel shape size
    TimeLevels<LADY_MESH_INDEX*, 2> idx;    //>! parcel mesh index
public:
    Parcel(int numDim);
    virtual ~Parcel();

    Parcel& operator=(const Parcel &other);
    
    /**
     *  Set the ID of tracer.
     *
     *  @param ID the given ID.
     */
    virtual void setID(int ID);

    /**
     *  Get the ID of tracer.
     *
     *  @return The ID.
     */
    virtual int getID() const;

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

    /**
     *  Get the shape function value for a given body coordinate.
     *
     *  @param timeIdx the time level index.
     *  @param y       the body coordinate.
     *
     *  @return The shape function value.
     */
    double getShapeFunction(const TimeLevelIndex<2> &timeIdx,
                            const LADY_BODY_COORD &y);

    void updateShapeSize(const LADY_DOMAIN &domain,
                         const TimeLevelIndex<2> &timeIdx);

    const vec::fixed<2>& getShapeSize(const TimeLevelIndex<2> &timeIdx) const;
};

}

#endif