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
    Foam::immobile
    
Description
    Transport model for stationary species (e.g., heavy ions in fast 
    timescale simulations). 
    
    The density changes only via time-evolution and source terms:
    dn/dt = Sources

    The model is activated by specifying:

        transportModel immobile;

    in the corresponding species dictionary.

SourceFiles
    immobile.C

\*---------------------------------------------------------------------------*/

#ifndef immobile_H
#define immobile_H

#include "plasmaTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class immobile Declaration
\*---------------------------------------------------------------------------*/

class immobile
:
    public plasmaTransportModel
{
private:

    //- Disallow default bitwise copy construct
    immobile(const immobile&) = delete;
    void operator=(const immobile&) = delete;
    
public:

    //- Runtime type information
    TypeName("immobile");

    // Constructors

        //- Construct from components
        immobile
        (
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            plasmaSpecies& species,
            const label specieIndex
        );

    //- Destructor
    virtual ~immobile() = default;

    // Member Functions
        
        //- No tranport coefficients to update
        virtual void correct() override;

        //- Return the dn/dt = 0 equation
        virtual tmp<fvScalarMatrix> nEqn() const override;

        //- No fluxes for immobile species — FatalError if called
        virtual void updateFluxes
        (
            const fvScalarMatrix& nEqnMatrix,
            surfaceScalarField& convectiveFlux,
            surfaceScalarField& diffusiveFlux,
            surfaceScalarField& particleFlux
        ) const override;

        //- Zero conductivity
        virtual tmp<volScalarField> electricalConductivity() const override;

        //- Zero diffusive source
        virtual tmp<volScalarField> diffusiveChargeSource() const override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
