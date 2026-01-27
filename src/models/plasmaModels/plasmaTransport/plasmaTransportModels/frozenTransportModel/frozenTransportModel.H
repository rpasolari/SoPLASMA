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
    Foam::frozenTransportModel
    
Description
    Transport model for stationary species (e.g., heavy ions in fast 
    timescale simulations). 
    
    The density changes only via time-evolution and source terms:
    dn/dt = Sources

    The model is activated by specifying:

        transportModel frozen;

    in the corresponding species dictionary.

SourceFiles
    frozenTransportModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_frozenTransportModel_H
#define Foam_frozenTransportModel_H

#include "plasmaTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class frozenTransportModel Declaration
\*---------------------------------------------------------------------------*/

class frozenTransportModel
:
    public plasmaTransportModel
{
private:

public:

    //- Runtime type information
    TypeName("frozen");

    // Constructors

        //- Construct from components
        frozenTransportModel
        (
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            const plasmaSpecies& species,
            const label specieIndex,
            const volVectorField& E,
            const volScalarField& ePotential
        );

    // Member Functions
        
        //- Correct model
        virtual void correct() override;

        //- Return the transport equation (Matrix)
        virtual tmp<fvScalarMatrix> nEqn() const override;

        //- Return the total particle flux (Advection + Diffusion) (1/s)
        virtual tmp<surfaceScalarField> particleFlux() const override;

        //- Return drift velocity
        virtual const volVectorField& driftVelocity() const override;
        
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
