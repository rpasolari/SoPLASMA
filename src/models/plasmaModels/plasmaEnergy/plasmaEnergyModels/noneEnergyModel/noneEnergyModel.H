/*---------------------------------------------------------------------------*\
License
    This file is part of the foamPlasmaToolkit.

    The foamPlasmaToolkit is not part of OpenFOAM but is developed using the
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
    Foam::noneEnergyModel
    
Description
    Energy model for species where energy tracking is disabled.
    
    This model is a placeholder that performs no calculations and owns 
    no temperature fields. It is typically used for background components 
    or species whose energy state is irrelevant to the simulation.

    The model is activated by specifying:

        energyModel none;

    in the corresponding species dictionary.

SourceFiles
    noneEnergyModel.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_noneEnergyModel_H
#define Foam_noneEnergyModel_H

#include "plasmaEnergyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class noneEnergyModel Declaration
\*---------------------------------------------------------------------------*/

class noneEnergyModel
:
    public plasmaEnergyModel
{
private:

public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from components
        noneEnergyModel
        (
            const word& modelName,
            const dictionary& dict,
            const fvMesh& mesh,
            const plasmaSpecies& species,
            const label specieIndex,
            const volVectorField& E
        );

    // Member Functions

        //- The correct does nothing
        virtual void correct() override;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
