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
    Foam::plasmaEnergy

Description
    High-level manager for plasma energy models.

    This class owns and manages one energy model per plasma species (mainly 
    for electrons). It does not implement any energy physics itself. Instead it
    serves as a controller that:

        - Retrieves species-level energy configuration
        - Constructs the appropriate energy model (LFA, LEA, etc.)
        - Provides temperature fields to the transport models

    The plasmaEnergy object is typically created once in the solver and used
    during the solution loop to update the energy behavior of the plasma
    species.

SourceFiles
    plasmaEnergyI.H
    plasmaEnergy.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_plasmaEnergy_H
#define Foam_plasmaEnergy_H

#include "fvCFD.H"
#include "fvMesh.H"
#include "PtrList.H"

#include "plasmaSpecies.H"
#include "plasmaEnergyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class plasmaEnergy Declaration
\*---------------------------------------------------------------------------*/

class plasmaEnergy
:
    public regIOobject
{
private:
    //- Flag to solve background gas energy equation
    Switch solveGasEnergy_;

    //- Flag to track if gas temperature is a field (solved or fixed)
    Switch isGasTempField_;

    //- The scalar value if isGasTempField_ == false
    dimensionedScalar TgasValue_;

    //- The volScalarField if isGasTempField_ == true
    autoPtr<volScalarField> TgasFieldPtr_;

    //- Reference to mesh
    const fvMesh& mesh_;

    //- Reference to species database
    plasmaSpecies& species_;

    //- Reference to the electric field
    const volVectorField& E_;

    //- Per-species energy models
    PtrList<plasmaEnergyModel> energyModels_;

    //- A gas temperature field for shared across the 

    //- Disallow copying
    plasmaEnergy(const plasmaEnergy&) = delete;
    void operator=(const plasmaEnergy&) = delete;

    //- Construct per-species models
    void constructModels();

public:

    //- Runtime type information
    TypeName("plasmaEnergy");

    // Constructors

        //- Construct from plasmaSpecies and mesh
        plasmaEnergy
        (
            plasmaSpecies& species,
            const fvMesh& mesh,
            const volVectorField& E
        );

    //- Destructor
    virtual ~plasmaEnergy() {}

    // Member functions

        //- Correct all models
        void correct();

        // //- Return total total Joule heating [W/m^3]
        // tmp<volScalarField> totalJouleHeating() const;

        virtual bool writeData(Ostream& os) const;

    // Access

        //- Mesh and species metadata access
        inline const fvMesh& mesh() const { return mesh_; }
        inline const plasmaSpecies& species() const { return species_; }

        //- Model access
inline plasmaEnergyModel& model(const label i)
        {
            return energyModels_[i];
        }

        inline const plasmaEnergyModel& model(const label i) const
        {
            return energyModels_[i];
        }

        //- Temperature field access for species i
inline tmp<volScalarField> T(const label i) const
        {
            return energyModels_[i].T();
        }

        //- Return the background gas temperature as a tmp field
        tmp<volScalarField> Tgas() const;

        //- Return the background gas temperature as a dimensionedScalar
        const dimensionedScalar& TgasValue() const;

//- Return true if gas temperature is a single uniform value (not a field)
inline bool isGasUniform() const { return !isGasTempField_; }

//- Return true if gas temperature is NOT being solved (remains constant in time)
inline bool isGasConstant() const { return !solveGasEnergy_; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "plasmaEnergyI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
