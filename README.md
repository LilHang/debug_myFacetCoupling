A user-defined `evalSourcesFromOutside` function.

 I just put it  directly in  `.../facet/box/couplingmanager.hh  `(dumux source code)  for box scheme (below `evalSourcesFromBulk`).

```c++
// header files copy from `../box/darcylaw.hh` in dumux source code
#include <vector>
#include <cmath>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/float_cmp.hh>

#include <dumux/common/parameters.hh>

#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/flux/cvfe/darcyslaw.hh>


NumEqVector<lowDimId> evalSourcesFromOutside(const Element<lowDimId>& element,
                                            const FVElementGeometry<lowDimId>& fvGeometry,
                                            const ElementVolumeVariables<lowDimId>& elemVolVars,
                                            const SubControlVolume<lowDimId>& scv,
                                            const double& transferCoeff)
{
    // make sure the this is called for the element of the context
    assert(this->problem(lowDimId).gridGeometry().elementMapper().index(element) == lowDimContext_.elementIdx);

    NumEqVector<lowDimId> sources(0.0);

    // find current coupled {LowDimIndexType:LowDimCouplingData}
    const auto& map = couplingMapperPtr_->couplingMap(lowDimGridId, bulkGridId);
    auto it = map.find(lowDimContext_.elementIdx);
    if (it == map.end())
        return sources;

    assert(lowDimContext_.isSet);

    // loop over all embedments (coupled bulk scvs)
    for (unsigned int i = 0; i < it->second.embedments.size(); ++i)
    {
        const auto& embedment = it->second.embedments[i];

        // list of scvfs in the bulk domain whose fluxes enter this scv
        // if low dim domain uses tpfa, this is all scvfs lying on this element
        // if it uses box, it is the one scvf coinciding with the given scv
        const auto& coincidingScvfs = embedment.second;
        const auto& scvfList = lowDimUsesBox ? std::vector<GridIndexType<lowDimId>>{ coincidingScvfs[scv.localDofIndex()] }
                                                : coincidingScvfs;

        // loop over all scvfs in the list
        NumEqVector<bulkId> coupledFluxes(0.0);
        for (const auto& scvfIdx : scvfList)
        {
            // get facetVolVars
            const auto& facetVolVars = this->getLowDimVolVars(this->problem(bulkId).gridGeometry().element(embedment.first),
                                                                (*(lowDimContext_.bulkFvGeometries[i])).scvf(scvfIdx));
            // compute outflow flux
            if (facetVolVars.pressure(0) > 1e5)
            {
                auto scvf = (*(lowDimContext_.bulkFvGeometries[i])).scvf(scvfIdx);
                const auto a = facetVolVars.extrusionFactor();
                auto gradP = scvf.unitOuterNormal();
                gradP *= 0.5 * a;
                gradP /= gradP.two_norm2();
                gradP *= (1e5 - facetVolVars.pressure(0));

                // 1000:rho, 1e-3:viscosity
                double conductance = 1000*facetVolVars.permeability()*transferCoeff/1e-3;
                // first 1.0 refer to positive sign, second 1.0 refer to the extrusion factor of the scvf
                coupledFluxes[0] += 1.0 * scvf.area() * 1.0
                                    * vtmv(scvf.unitOuterNormal(), conductance, gradP);
            }
        }

        sources += coupledFluxes;
    }

    return sources;
}

```

