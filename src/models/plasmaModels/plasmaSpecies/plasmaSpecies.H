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
    Foam::plasmaSpecies
    
Description
    Core class for plasma species data management.

    This class provides the minimal data structures required to represent
    plasma species in OpenFOAM-based plasma simulations. It defines which
    species exist, their basic properties (e.g. charge number), and
    constructs and owns the corresponding volume fields such as number
    density.

    plasmaSpecies does not implement any physical models (transport,
    reactions, electron energy, etc.). Instead, it acts as a central storage
    and access point for all species-related quantities. Other physics
    modules operate on these fields but do not create their own copies.

    Only a single plasmaSpecies object should exist on a given mesh.
    It is responsible for creating the species fields; additional modules
    simply reference these fields to ensure consistency across the solver.

    This class forms the foundation upon which transport models, reaction
    models, and electrostatic models build their physics algorithms.

SourceFiles
    plasmaSpeciesI.H
    plasmaSpecies.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_plasmaSpecies_H
#define Foam_plasmaSpecies_H

#include "IOdictionary.H"
#include "fvMesh.H"
#include "volFields.H"
#include "HashTable.H"

#include "foamPlasmaToolkitConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class basicPlasmaSpecies Declaration
\*---------------------------------------------------------------------------*/

class plasmaSpecies
:
    public IOdictionary
{
private:

    //- Disallow copy
    plasmaSpecies(const plasmaSpecies&) = delete;
    void operator=(const plasmaSpecies&) = delete;

protected:

    // Protected Data

        //- Mesh reference
        const fvMesh& mesh_;

        //- Active species names
        wordList speciesNames_;

        //- Active species index-name mapping
        HashTable<label> speciesID_;

        //- Number of active species
        label nSpecies_;

        //- Charge number Z_i per species (dimensionless:: -1, +1, +2, etc)
        List<scalar> speciesChargeNumber_;

        //- Charge per species
        PtrList<dimensionedScalar> speciesCharge_;

        //- Mass per species
        PtrList<dimensionedScalar> speciesMass_;

        //- Number density fiels n_i [1/m3] for each active species
        PtrList<volScalarField> numberDensity_;

        //- Per-species property dictionaries
        HashTable<dictionary> speciesDicts_;

        //- Default species properties
        dictionary defaultDict_;

        //- Background gas dict
        dictionary backgroundGasDict_;

        //- Label for electron specie
        label electronSpeciesID_;

        //- Label list for ions (positive + negative)
        labelList ionSpeciesIDs_;

        //- Label list for positive ions only
        labelList positiveIonSpeciesIDs_;

        //- Label list for negative ions only
        labelList negativeIonSpeciesIDs_;

        //- Label list for active neutral species (e.g., excited species)
        labelList neutralSpeciesIDs_;

        //- Label list for charged species (ions + electrons)
        labelList chargedSpeciesIDs_;

        //- Label list for mobile species
        labelList mobileSpeciesIDs_;

        //- Label list for immobile species
        labelList immobileSpeciesIDs_;

        //- Species with a constant temperature field
        labelList constantTemperatureSpeciesIDs_;

        //- Species with a dynamic temperature field
        labelList dynamicTemperatureSpeciesIDs_;

        //- Species solving the energy equation
        labelList solveEnergyEquationSpeciesIDs_;

        //- Species with a field temperature (not just a dimensionedScalar)
        labelList fieldTemperatureSpeciesIDs_;

public:

    //- Runtime type information
    TypeName("plasmaSpecies");

    // Constructors

        //- Construct from mesh
        plasmaSpecies
        (
            const fvMesh& mesh
        );

    //- Destructor
    virtual ~plasmaSpecies() {}

    // Access

        inline const fvMesh& mesh() const;
        inline const wordList& speciesNames() const;
        inline label nSpecies() const;

        inline const List<scalar>& speciesChargeNumber() const;
        inline scalar speciesChargeNumber(label i) const;
        inline const PtrList<dimensionedScalar>& speciesCharge() const;
        inline const dimensionedScalar& speciesCharge(label i) const;
        inline const PtrList<dimensionedScalar>& speciesMass() const;
        inline const dimensionedScalar& speciesMass(label i) const;

        //- Get species index from name
        inline label speciesID(const word& name) const;

        //- Number density access (single specie and whole list)
        inline const volScalarField& numberDensity(label i) const;
        inline volScalarField& numberDensity(label i);
        inline const PtrList<volScalarField>& numberDensity() const;
        inline PtrList<volScalarField>& numberDensity();

        //- Background gas dictionary access
        inline const dictionary& backgroundGasDict() const;

        //- Species dictionaries access
        inline const dictionary& speciesDict(label i) const;
        inline const dictionary& speciesDict(const word& name) const; 
        inline const dictionary& defaultDict() const; 

        //- Group Access
        inline label electronSpeciesID() const;
        inline const labelList& ionSpeciesIDs() const;
        inline const labelList& positiveIonSpeciesIDs() const;
        inline const labelList& negativeIonSpeciesIDs() const;
        inline const labelList& neutralSpeciesIDs() const;
        inline const labelList& chargedSpeciesIDs() const;
        inline const labelList& mobileSpeciesIDs() const;
        inline const labelList& immobileSpeciesIDs() const;
        inline const labelList& constantTemperatureSpeciesIDs() const;
        inline const labelList& dynamicTemperatureSpeciesIDs() const;
        inline const labelList& solveEnergyEquationSpeciesIDs() const;
        inline const labelList& fieldTemperatureSpeciesIDs() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "plasmaSpeciesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
