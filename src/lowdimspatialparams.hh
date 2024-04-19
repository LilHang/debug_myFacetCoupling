// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetTests
 * \brief The spatial parameters for the richards facet coupling test.
 */
#ifndef DUMUX_TEST_LOWDIM_RICHARDS_SPATIALPARAMS_HH
#define DUMUX_TEST_LOWDIM_RICHARDS_SPATIALPARAMS_HH

#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup FacetTests
 * \brief The spatial parameters for the richards facet coupling test.
 */
template< class GridGeometry, class Scalar >
class LowDimSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP< GridGeometry, Scalar, LowDimSpatialParams<GridGeometry, Scalar> >
{
    using ThisType = LowDimSpatialParams< GridGeometry, Scalar >;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP< GridGeometry, Scalar, ThisType >;

    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    // lihang add start
    using Grid = typename GridView::Grid;
    // lihang add end
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;


public:
    //! Export the type used for permeabilities
    using PermeabilityType = Scalar;

    LowDimSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                        // lihang add start
                        std::shared_ptr<const Dumux::GridData<Grid>> gridData,
                        // lihang add end
                        const std::string& paramGroup)
    : ParentType(gridGeometry)
    // lihang add start
    , gridDataPtr_(gridData)
    // lihang add end
    {
        permeability_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Permeability");
        aperture_ = getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture");
    }

    //! Function for defining the (intrinsic) permeability \f$[m^2]\f$.
    template< class ElementSolution >
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        return permeability_;
    }

     //! Returns the porosity
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }



    /*!
     * \brief Returns the temperature at the domain at the given position
     * \param globalPos The position in global coordinates where the temperature should be specified
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        return 283.15;
    }

    //! Set the aperture as extrusion factor.
    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
    {
        // We treat the fractures as lower-dimensional in the grid,
        // but we have to give it the aperture as extrusion factor
        // such that the dimensions are correct in the end.
        return aperture_;
    }

    // lihang add start
    // returns the domain marker for an element
    int getElementDomainMarker(const Element& element) const
    {
        return gridDataPtr_->getElementDomainMarker(element);
    }
    // lihang add end


private:
    // lihang add start
    // pointer to the grid data (contains domain markers)
    std::shared_ptr<const Dumux::GridData<Grid>> gridDataPtr_;
    // lihang add end

    PermeabilityType permeability_;
    Scalar aperture_;
};

} // end namespace Dumux

#endif
