// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetTests
 * \brief The properties for the bulk domain in the rinchards facet coupling test.
 */
#ifndef DUMUX_TEST_FACETCOUPLING_RICHARDS_MATRIX_PROPERTIES_HH
#define DUMUX_TEST_FACETCOUPLING_RICHARDS_MATRIX_PROPERTIES_HH

#include <dumux/material/components/constant.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/cctpfa.hh>

#include <dumux/multidomain/facet/box/properties.hh>
#include <dumux/multidomain/facet/cellcentered/tpfa/properties.hh>
#include <dumux/multidomain/facet/cellcentered/mpfa/properties.hh>
#include <dumux/multidomain/facet/couplingmapper.hh>
#include <dumux/multidomain/facet/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

// Matrix sub-problem
// the spatial parameters (permeabilities, material parameters etc.)
#include "matrixspatialparams.hh"
// we use alu grid for the discretization of the matrix domain
#include <dune/alugrid/grid.hh>

// LowDim sub-problem
// the spatial parameters (permeabilities, material parameters etc.)
#include "lowdimspatialparams.hh"
// we use foam grid for the discretization of the fracture domain
// as this grid manager is able to represent network/surface grids
#include <dune/foamgrid/foamgrid.hh>


#include "problem_matrix.hh"
#include "problem_lowdim.hh"

// default for the bulk grid type
#ifndef MATRIXGRIDTYPE
#define MATRIXGRIDTYPE Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>
#endif

#ifndef LOWDIMGRIDTYPE
#define LOWDIMGRIDTYPE Dune::FoamGrid<1, 2>
#endif

namespace Dumux::Properties {

// create the type tag nodes
namespace TTag {
struct OnePMatrix { using InheritsFrom = std::tuple<OneP>; };
struct OnePMatrixTpfa { using InheritsFrom = std::tuple<CCTpfaFacetCouplingModel, OnePMatrix>; };
// struct OnePMatrixMpfa { using InheritsFrom = std::tuple<CCMpfaFacetCouplingModel, OnePMatrix>; };
struct OnePMatrixBox { using InheritsFrom = std::tuple<BoxFacetCouplingModel, OnePMatrix>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePMatrix> { using type = MATRIXGRIDTYPE; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePMatrix> { using type = OnePMatrixProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePMatrix>
{
    using type = MatrixSpatialParams< GetPropType<TypeTag, Properties::GridGeometry>,
                                    GetPropType<TypeTag, Properties::Scalar> >;
};

// the fluid system
    template<class TypeTag>
    struct FluidSystem<TypeTag, TTag::OnePMatrix>
    {
    private:
        using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    public:
        using type = FluidSystems::OnePLiquid< Scalar, Components::Constant<1, Scalar> >;
    };


// create the type tag nodes
namespace TTag {
struct OnePLowDim { using InheritsFrom = std::tuple<OneP>; };
struct OnePLowDimTpfa { using InheritsFrom = std::tuple<OnePLowDim, CCTpfaModel>; };

// we need an additional type tag for the test using box in the bulk domain
//struct OnePLowDimMpfa { using InheritsFrom = std::tuple<OnePLowDim, CCTpfaModel>; };
struct OnePLowDimBox { using InheritsFrom = std::tuple<OnePLowDim, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePLowDim> { using type = LOWDIMGRIDTYPE; };
// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePLowDim> { using type = OnePLowDimProblem<TypeTag>; };
// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePLowDim>
{
    using type = LowDimSpatialParams< GetPropType<TypeTag, Properties::GridGeometry>,
                                    GetPropType<TypeTag, Properties::Scalar> >;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePLowDim>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid< Scalar, Components::Constant<1, Scalar> >;
};


// obtain/define some types to be used below in the property definitions and in main
template< class MatrixTypeTag, class LowDimTypeTag >
class TestTraits
{
    using MatrixFVGridGeometry = GetPropType<MatrixTypeTag, Properties::GridGeometry>;
    using LowDimFVGridGeometry = GetPropType<LowDimTypeTag, Properties::GridGeometry>;
public:
    using MDTraits = Dumux::MultiDomainTraits<MatrixTypeTag, LowDimTypeTag>;
    using CouplingMapper = Dumux::FacetCouplingMapper<MatrixFVGridGeometry, LowDimFVGridGeometry>;
    using CouplingManager = Dumux::FacetCouplingManager<MDTraits, CouplingMapper>;
};

// set cm property in the sub-problems
using TpfaTraits = TestTraits<TTag::OnePMatrixTpfa, TTag::OnePLowDimTpfa>;
using BoxTraits = TestTraits<TTag::OnePMatrixBox, TTag::OnePLowDimBox>;
template<class TypeTag> struct CouplingManager<TypeTag, TTag::OnePMatrixTpfa> { using type = typename TpfaTraits::CouplingManager; };
template<class TypeTag> struct CouplingManager<TypeTag, TTag::OnePLowDimTpfa> { using type = typename TpfaTraits::CouplingManager; };
template<class TypeTag> struct CouplingManager<TypeTag, TTag::OnePMatrixBox> { using type = typename BoxTraits::CouplingManager; };
template<class TypeTag> struct CouplingManager<TypeTag, TTag::OnePLowDimBox> { using type = typename BoxTraits::CouplingManager; };

} // end namespace Dumux::Properties

#endif
