/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2026
        Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::driftDiffusion
    
Description
    Drift–diffusion transport model for plasma species.

    This model implements the standard drift–diffusion formulation, combining
    mobility-driven drift in the electric field with concentration-driven
    diffusion. It updates the species fluxes and densities according to the
    provided mobility and diffusion coefficients.

    The model is activated by specifying:

        transportModel driftDiffusion;

    in the corresponding species dictionary.

SourceFiles
    driftDiffusion.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_driftDiffusion_H
#define Foam_driftDiffusion_H

#include "plasmaTransportModel.H"
#include "plasmaMobilityModel.H"
#include "plasmaDiffusivityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class driftDiffusion Declaration
\*---------------------------------------------------------------------------*/

class driftDiffusion
:
    public plasmaTransportModel
{
private:

    //- Disallow default bitwise copy construct
    driftDiffusion(const driftDiffusion&) = delete;
    //- Disallow default bitwise assignment
    void operator=(const driftDiffusion&) = delete;

    //- Mobility/Diffusivity submodels
    autoPtr<plasmaMobilityModel> mobilityModel_;
    autoPtr<plasmaDiffusivityModel> diffusivityModel_;
    
    //- Configuration settings
    word fluxScheme_;

    // Cached fluxes

        // Cached convective flux surface field (phi = Z*mu*phiE)
        mutable autoPtr<surfaceScalarField> cachedPhi_;
    
    // Fluxes

        //- Convective particle flux [1/s]
        mutable surfaceScalarField convectiveFlux_;

        //- Diffusive particle flux [1/s]
        mutable surfaceScalarField diffusiveFlux_;

        //- Total particle flux = convective + diffusive [1/s]
        mutable surfaceScalarField particleFlux_;

    //- Build convective flux
    tmp<surfaceScalarField> convectivePhi() const;

    //- Construct mobility and diffusivity models
    void constructModels();

    //- Allocate the flux fields
    void initFluxFields();

    //- Split SG flux into convective and diffusive components
    void splitSGFlux
    (
        const surfaceScalarField& phi,
        const volScalarField& D,
        const volScalarField& n,
        surfaceScalarField& convFlux,
        surfaceScalarField& diffFlux
    ) const;

public:

    //- Runtime type information
    TypeName("driftDiffusion");

    // Constructors

        //- Construct from components
        driftDiffusion
        (
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            plasmaSpecies& species,
            const label specieIndex
        );

    //- Destructor

    virtual ~driftDiffusion() = default;

    // Member Functions
        
        //- Update mu and D fields from evaluators
        virtual void correct() override;

        //- Return number density equation matrix
        virtual tmp<fvScalarMatrix> nEqn() const override;

        //- Update convective, diffusive, and total particle fluxes
        virtual void updateFluxes
        (
            const fvScalarMatrix& nEqnMatrix,
            surfaceScalarField& convectiveFlux,
            surfaceScalarField& diffusiveFlux,
            surfaceScalarField& particleFlux
        ) const override;

        //- Return electrical conductivity [S/m]
        virtual tmp<volScalarField> electricalConductivity() const override;

        //- Return diffusive charge source [C/m3/s]
        virtual tmp<volScalarField> diffusiveChargeSource() const override;

    // Access

        //- Return the convective flux surface field
        const surfaceScalarField& convectiveFlux() const
        {
            return convectiveFlux_;
        }

        //- Return the diffusive flux surface field
        const surfaceScalarField& diffusiveFlux() const
        {
            return diffusiveFlux_;
        }

        //- Return the total flux surface field
        const surfaceScalarField& particleFlux() const
        {
            return particleFlux_;
        }

        //- Return the mobility model
        const plasmaMobilityModel& mobility() const
        {
            return *mobilityModel_;
        }

        //- Return the diffusivity model
        const plasmaDiffusivityModel& diffusivity() const
        {
            return *diffusivityModel_;
        }

        //- Return the flux scheme
        const word& fluxScheme() const { return fluxScheme_; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
