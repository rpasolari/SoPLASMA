/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2025
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
    //- Mobility/Diffusivity submodels
    autoPtr<plasmaMobilityModel> mobilityModel_;
    autoPtr<plasmaDiffusivityModel> diffusivityModel_;
    
    //- Configuration settings
    word advection_;
    word fluxScheme_;
    bool isExplicit_;

    //- Construct mobility and diffusivity models
    void constructModels();

    //- Disallow default bitwise copy construct
    driftDiffusion(const driftDiffusion&) = delete;

    //- Disallow default bitwise assignment
    void operator=(const driftDiffusion&) = delete;

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
            const plasmaSpecies& species,
            const label specieIndex,
            const volVectorField& E,
            const surfaceScalarField& phiE
        );

    //- Destructor

    virtual ~driftDiffusion() = default;

    // Member Functions
        
        //- Return the transport equation (Matrix)
        virtual tmp<fvScalarMatrix> nEqn() const override;

        //- Update the particle flux in place (1/s)
        virtual void updateParticleFlux
        (
            surfaceScalarField& flux
        ) const override;

        //- Return electric conduction coefficient
        virtual tmp<volScalarField> electricalConductivity() const override;

        //- Return electric diffusion charge
        virtual tmp<volScalarField> diffusiveChargeSource() const override;

        //- Return the flux scheme
        const word& fluxScheme() const { return fluxScheme_; }

        const plasmaMobilityModel& mobility() const
        {
            return *mobilityModel_;
        }

        const plasmaDiffusivityModel& diffusivity() const
        {
            return *diffusivityModel_;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
