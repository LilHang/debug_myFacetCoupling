// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FacetTests
 * \brief The problem for the bulk domain in the richards facet coupling test.
 */
#ifndef DUMUX_TEST_FACETCOUPLING_RICHARDS_BULKPROBLEM_HH
#define DUMUX_TEST_FACETCOUPLING_RICHARDS_BULKPROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>

// lihang add start
#include <vector>
#include <cmath>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/common/parameters.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/flux/cvfe/darcyslaw.hh>
// lihang add end

namespace Dumux {

/*!
 * \ingroup FacetTests
 * \brief Test problem for the incompressible richards model
 *        with coupling across the bulk grid facets.
 */
template<class TypeTag>
class OnePMatrixProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename GridVariables::Scalar;

    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    // some indices
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    // static constexpr int dimWorld = GridView::dimensionworld;

    // lihang add start
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    // lihang add end

public:
    OnePMatrixProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                    std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                    std::shared_ptr<CouplingManager> couplingManager,
                    const std::string& paramGroup = "")
    : ParentType(gridGeometry, spatialParams, paramGroup)
    , couplingManagerPtr_(couplingManager)
    {
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" +
                         getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");

        // lihang add start
        // read parameters from input file
        transferCoeff_ = getParamFromGroup<Scalar>("LowDim", "SpatialParams.TransferCoefficient", 1.0);
        const auto boundaryHeadInCm = getParamFromGroup<Scalar>("Matrix", "BoundaryHeadInCm");
        boundaryPressure_ = 1.0e5 + 1000.0*9.81*boundaryHeadInCm*0.01;
        // lihang add end
    }

    //! The problem name.
    const std::string& name() const
    { return problemName_; }

    //! Specifies the type of boundary condition at a given position.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes bcTypes;
        if (globalPos[0] < 0.215 + 1e-6)
            bcTypes.setAllNeumann();
        else
          bcTypes.setAllDirichlet();
        return bcTypes;
    }

    //! Specifies the type of interior boundary condition at a given position.
    BoundaryTypes interiorBoundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann();
        return bcTypes;
    }

    //! Evaluates the source term at a given position.
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

    //! Evaluates the Dirichlet boundary condition for a given position.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables  values(0.0);
        values[Indices::pressureIdx] = boundaryPressure_;
//        std::cout << "dirichletAtPos: " << globalPos << std::endl;
        return values;
    }

    //! Evaluates the Neumann boundary condition for a boundary segment.
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

    //! Evaluates the initial conditions.
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables  values(0.0);
        values[Indices::pressureIdx] = 1.0e5;
        return values;
    }

    // reference pressure
    Scalar nonwettingReferencePressure() const
    { return 1.0e5; };

    //! Returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

    // lihang add start
    // for output purpose only
    template<class SolutionVector, class GridVariables, class Assembler, class MatrixDomainId>
    Scalar computeOutFlow(const SolutionVector& sol, const GridVariables& gridVars, const MatrixDomainId& matrixDomainId,
                          const Assembler& assembler, bool verbose = true) const
    {
        Scalar rate = 0.0;

        auto fvGeometry = localView(this->gridGeometry());


        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
          couplingManagerPtr_->bindCouplingContext(matrixDomainId, element, assembler);
            fvGeometry.bindElement(element);

            for (const auto& scvf : scvfs(fvGeometry))
                if (couplingManagerPtr_->isOnInteriorBoundary(element, scvf))
                {
                    const auto& facetVolVars = couplingManagerPtr_->getLowDimVolVars(element, scvf);
                    if (facetVolVars.pressure(0) > 1e5)
                    {
                        const auto a = facetVolVars.extrusionFactor();
                        auto gradP = scvf.unitOuterNormal();
                        gradP *= 0.5 * a;
                        gradP /= gradP.two_norm2();
                        gradP *= (1e5 - facetVolVars.pressure(0));

                        // 1000:rho, 1e-3:viscosity
                        double conductance = 1000*facetVolVars.permeability()*transferCoeff_/1e-3;
                        rate += 1.0 * scvf.area() * 1.0
                                            * vtmv(scvf.unitOuterNormal(), conductance, gradP);
                    }
                }
        }



        if (verbose)
            std::cout << Fmt::format("Actual rate: {:.5g} (m^2/s)\n", rate);

        return rate;
    }

    // for output purpose only
    template<class SolutionVector, class GridVariables, class Assembler, class MatrixDomainId>
    Scalar computeInterfaceFlux(const SolutionVector& sol, const GridVariables& gridVars, const MatrixDomainId& matrixDomainId,
                          const Assembler& assembler, bool verbose = true) const
    {
        Scalar rate = 0.0;

        auto fvGeometry = localView(this->gridGeometry());
        auto elemVolVars = localView(gridVars.curGridVolVars());
        auto elemFluxVarsCache = localView(gridVars.gridFluxVarsCache());

        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            couplingManagerPtr_->bindCouplingContext(matrixDomainId, element, assembler);
            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, sol);
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

            for (const auto& scvf : scvfs(fvGeometry))
                if (couplingManagerPtr_->isOnInteriorBoundary(element, scvf))
                {
                    FluxVariables fluxVars;
                    fluxVars.init(*this, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

                    // 0 refer to the first phase
                    auto upwindTerm = [](const auto& volVars)
                    { return volVars.density(0)*volVars.mobility(0); };
                    Scalar flux = fluxVars.advectiveFlux(0, upwindTerm);
                    rate += flux;

                }
        }

        if (verbose)
            std::cout << Fmt::format("Actual rate: {:.5g} (m^2/s)\n", rate);

        return rate;
    }
    // lihang add end

private:
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    std::string problemName_;
    Scalar transferCoeff_;
    Scalar boundaryPressure_;
};

} // end namespace Dumux

#endif
