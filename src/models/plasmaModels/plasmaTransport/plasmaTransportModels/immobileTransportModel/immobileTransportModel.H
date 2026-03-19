/*---------------------------------------------------------------------------*\
License
    This file is part of the foamPlasmaToolkit.

    The foamPlasmaToolkit is not part of OpenFOAM but is developed using the
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
    Foam::immobileTransportModel
    
Description
    Transport model for stationary species (e.g., heavy ions in fast 
    timescale simulations). 
    
    The density changes only via time-evolution and source terms:
    dn/dt = Sources

    The model is activated by specifying:

        transportModel immobile;

    in the corresponding species dictionary.

SourceFiles
    immobileTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_immobileTransportModel_H
#define Foam_immobileTransportModel_H

#include "plasmaTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class immobileTransportModel Declaration
\*---------------------------------------------------------------------------*/

class immobileTransportModel
:
    public plasmaTransportModel
{
public:

    //- Runtime type information
    TypeName("immobile");

    // Constructors

        //- Construct from components
        immobileTransportModel
        (
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            const plasmaSpecies& species,
            const label specieIndex,
            const volVectorField& E,
            const surfaceScalarField& phiE
        );

    // Member Functions
        
        //- Return the transport equation (Matrix)
        virtual tmp<fvScalarMatrix> nEqn() const override;

        //- Update particle flux
        virtual void updateParticleFlux(surfaceScalarField& flux) const override;

        //- Return electric conduction coefficient
        virtual tmp<volScalarField> electricalConductivity() const override;

        //- Return electric diffusion charge
        virtual tmp<volScalarField> diffusiveChargeSource() const override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
