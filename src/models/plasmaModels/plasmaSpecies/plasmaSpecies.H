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
#include "PtrList.H"
#include "labelList.H"
#include "wordList.H"

#include "foamPlasmaToolkitConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class plasmaSpecies Declaration
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

        //- Active spcecies

            //- Number of active species
            label nSpecies_;

            //- Active species names
            wordList speciesNames_;

            //- Active species index-name mapping
            HashTable<label> speciesIDs_;

            //- Mass per active species [kg]
            PtrList<dimensionedScalar> speciesMasses_;

            //- Charge per active species [C]
            PtrList<dimensionedScalar> speciesCharges_;

            //- Charge number Z_i per active species (-1, +1, etc)
            List<scalar> speciesChargeNumbers_;

            //- Number density per active species [1/m3]
            PtrList<volScalarField> numberDensities_;

            //- Total charge density [C/m3]
            volScalarField chargeDensity_;

            //- Property dictionaries per active species
            HashTable<dictionary> speciesDicts_;

            //- Default active species properties
            dictionary defaultSpeciesDict_;

        //- Background gas

            //- Background gas name
            word backgroundName_;

            //- Background gas number density N [1/m3]
            dimensionedScalar backgroundDensity_;

            //- Total neutral number density field [1/m3] 
            //- Sum of nBackground_ + all active neutrals
            volScalarField totalNeutralDensity_;

            //- Background gas dictionary
            dictionary backgroundDict_;

        //- Groups

            //- Index for the electron species
            label electronSpeciesID_;

            //- Indices for all ions (positive + negative)
            labelList ionSpeciesIDs_;

            //- Indices for positive ions
            labelList positiveIonSpeciesIDs_;

            //- Indices for negative ions
            labelList negativeIonSpeciesIDs_;

            //- Indices for active neutrals
            labelList neutralSpeciesIDs_;

            //- Indices for all charged species (ions + electrons)
            labelList chargedSpeciesIDs_;

            //- Indices for mobile species
            labelList mobileSpeciesIDs_;

            //- Indices for immobile species
            labelList immobileSpeciesIDs_;

            //- Species with a constant (fixed) temperature value
            labelList constantTemperatureSpeciesIDs_;

            //- Species with a dynamic (non-fixed) temperature
            labelList dynamicTemperatureSpeciesIDs_;

            //- Species that follow the background gas temperature
            labelList followBackgroundTempSpeciesIDs_;

            //- Species for which we solve a dedicated energy equation
            labelList solveEnergySpeciesIDs_;

            //- Species that have a volScalarField for temperature
            labelList fieldTemperatureSpeciesIDs_;

            //- Species that participate in reactions
            labelList reactingSpeciesIDs_;

            //- Species that are inert (transport only, no chemical reactions)
            labelList nonReactingSpeciesIDs_;

public:

    //- Runtime type information
    TypeName("plasmaSpecies");

    // Constructors

        //- Construct from mesh
        explicit plasmaSpecies(const fvMesh& mesh);

    //- Destructor
    virtual ~plasmaSpecies() = default;

    // Member functions
        
        void update();

    // Access

        //- Metadata access
        inline const fvMesh& mesh() const;
        inline label nSpecies() const;
        inline const wordList& speciesNames() const;
        inline const word& speciesName(label i) const;

        //- Get species index from name
        inline label speciesID(const word& name) const;

        //- Physical properties access
        inline const List<scalar>& speciesChargeNumbers() const;
        inline scalar speciesChargeNumber(label i) const;
        inline const PtrList<dimensionedScalar>& speciesCharges() const;
        inline const dimensionedScalar& speciesCharge(label i) const;
        inline const PtrList<dimensionedScalar>& speciesMasses() const;
        inline const dimensionedScalar& speciesMass(label i) const;
        inline const PtrList<volScalarField>& numberDensities() const;
        inline PtrList<volScalarField>& numberDensities();
        inline const volScalarField& numberDensity(label i) const;
        inline volScalarField& numberDensity(label i);
        inline const volScalarField& chargeDensity() const;

        //- Dictionary access
        inline const dictionary& speciesDict(label i) const;
        inline const dictionary& speciesDict(const word& name) const; 
        inline const HashTable<dictionary>& speciesDicts() const;
        inline const dictionary& defaultSpeciesDict() const;

        //- Background gas access
        inline const word& backgroundName() const;
        inline const dimensionedScalar& backgroundDensity() const;
        inline const volScalarField& totalNeutralDensity() const;
        inline const dictionary& backgroundDict() const;

        //- Group access
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
        inline const labelList& followBackgroundTempSpeciesIDs() const;
        inline const labelList& solveEnergySpeciesIDs() const;
        inline const labelList& fieldTemperatureSpeciesIDs() const;
        inline const labelList& reactingSpeciesIDs() const;
        inline const labelList& nonReactingSpeciesIDs() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "plasmaSpeciesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
