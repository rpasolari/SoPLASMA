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
    Foam::gasTemperatureEnergyModel
    
Description
    Energy model for species which have the same temperature as the background
    gas.
    
    The model is activated by specifying:

        energyModel gasTemperature;

    in the corresponding species dictionary.

SourceFiles
    gasTemperatureEnergyModel.C

\*---------------------------------------------------------------------------*/

#ifndef gasTemperatureEnergyModel_H
#define gasTemperatureEnergyModel_H

#include "plasmaEnergyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class gasTemperatureEnergyModel Declaration
\*---------------------------------------------------------------------------*/

class gasTemperatureEnergyModel
:
    public plasmaEnergyModel
{
private:

    //- Disallow default bitwise copy construct
    gasTemperatureEnergyModel(const gasTemperatureEnergyModel&) = delete;

    //- Disallow default bitwise assignment
    void operator=(const gasTemperatureEnergyModel&) = delete;

public:

    //- Runtime type information
    TypeName("gasTemperature");

    // Constructors

        //- Construct from components
        gasTemperatureEnergyModel
        (
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            const plasmaSpecies& species,
            const label specieIndex,
            const volVectorField& E
        );

    // Member functions

        //- Correct model 
        virtual void correct() override;

        //- Return the energy equation (Matrix)
        virtual tmp<fvScalarMatrix> eEqn() const override;

        //- Returns true if the temperature is uniform (constant in space)
        virtual bool isUniform() const override;

        //- Returns true if the temperature is constant in time
        virtual bool isConstant() const override;

        //- Return a uniform field representation of the scalar T (K)
        virtual  tmp<volScalarField> T() const override;

        // Return the scalar value T (K)
        virtual const dimensionedScalar& Tvalue() const override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
