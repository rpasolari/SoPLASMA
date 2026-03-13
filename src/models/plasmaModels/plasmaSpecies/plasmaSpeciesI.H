/*---------------------------------------------------------------------------*\
  File: plasmaSpeciesI.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Inline member function definitions for the plasmaSpecies class.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef Foam_plasmaSpeciesI_H
#define Foam_plasmaSpeciesI_H

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Mesh and species metadata access
inline const Foam::fvMesh&
Foam::plasmaSpecies::mesh() const
{
    return mesh_;
}

inline const Foam::wordList&
Foam::plasmaSpecies::speciesNames() const
{
    return speciesNames_;
}

inline Foam::label
Foam::plasmaSpecies::nSpecies() const
{
    return nSpecies_;
}

inline const Foam::List<Foam::scalar>&
Foam::plasmaSpecies::speciesChargeNumber() const
{
    return speciesChargeNumber_;
}

inline Foam::scalar
Foam::plasmaSpecies::speciesChargeNumber(Foam::label i) const
{
    return speciesChargeNumber_[i];
}

inline const Foam::PtrList<Foam::dimensionedScalar>&
Foam::plasmaSpecies::speciesCharge() const
{
    return speciesCharge_;
}

inline const Foam::dimensionedScalar&
Foam::plasmaSpecies::speciesCharge(Foam::label i) const
{
    return speciesCharge_[i];
}

inline const Foam::PtrList<Foam::dimensionedScalar>&
Foam::plasmaSpecies::speciesMass() const
{
    return speciesMass_;
}

inline const Foam::dimensionedScalar&
Foam::plasmaSpecies::speciesMass(Foam::label i) const
{
    return speciesMass_[i];
}

inline Foam::label
Foam::plasmaSpecies::speciesID(const Foam::word& name) const
{
    return speciesID_[name];
}

//- Number density access
inline const Foam::volScalarField&
Foam::plasmaSpecies::numberDensity(Foam::label i) const
{
    return numberDensity_[i];
}

inline Foam::volScalarField&
Foam::plasmaSpecies::numberDensity(Foam::label i)
{
    return numberDensity_[i];
}

inline const Foam::PtrList<Foam::volScalarField>&
Foam::plasmaSpecies::numberDensity() const
{
    return numberDensity_;
}

inline Foam::PtrList<Foam::volScalarField>&
Foam::plasmaSpecies::numberDensity()
{
    return numberDensity_;
}


//- Species dictionaries access
inline const Foam::dictionary&
Foam::plasmaSpecies::speciesDict(Foam::label i) const
{
    return speciesDicts_[speciesNames_[i]];
}

inline const Foam::dictionary&
Foam::plasmaSpecies::speciesDict(const Foam::word& name) const
{
    return speciesDicts_[name];
}

inline const Foam::dictionary&
Foam::plasmaSpecies::defaultDict() const
{
    return defaultDict_;
}

//- Group Access
inline Foam::label
Foam::plasmaSpecies::electronSpeciesID() const
{
    return electronSpeciesID_;
}

inline const Foam::labelList&
Foam::plasmaSpecies::ionSpeciesIDs() const
{
    return ionSpeciesIDs_;
}

inline const Foam::labelList&
Foam::plasmaSpecies::positiveIonSpeciesIDs() const
{
    return positiveIonSpeciesIDs_;
}

inline const Foam::labelList&
Foam::plasmaSpecies::negativeIonSpeciesIDs() const
{
    return negativeIonSpeciesIDs_;
}

inline const Foam::labelList&
Foam::plasmaSpecies::neutralSpeciesIDs() const
{
    return neutralSpeciesIDs_;
}

inline const Foam::labelList&
Foam::plasmaSpecies::chargedSpeciesIDs() const
{
    return chargedSpeciesIDs_;
}

inline const Foam::labelList&
Foam::plasmaSpecies::mobileSpeciesIDs() const
{
    return mobileSpeciesIDs_;
}

inline const Foam::labelList&
Foam::plasmaSpecies::immobileSpeciesIDs() const
{
    return immobileSpeciesIDs_;
}

// inline const Foam::labelList&
// Foam::plasmaSpecies::nonIsothermalSpeciesIDs() const
// {
//     return nonIsothermalSpeciesIDs_;
// }
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
