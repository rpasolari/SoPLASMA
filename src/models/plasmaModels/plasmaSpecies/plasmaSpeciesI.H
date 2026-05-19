/*---------------------------------------------------------------------------*\
  File: plasmaSpeciesI.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Inline member function definitions for the plasmaSpecies class.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef Foam_plasmaSpeciesI_H
#define Foam_plasmaSpeciesI_H

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Metadata access
inline const Foam::fvMesh&
Foam::plasmaSpecies::mesh() const
{
    return mesh_;
}

inline Foam::label
Foam::plasmaSpecies::nSpecies() const
{
    return nSpecies_;
}

inline const Foam::wordList&
Foam::plasmaSpecies::speciesNames() const
{
    return speciesNames_;
}

inline const Foam::word&
Foam::plasmaSpecies::speciesName(Foam::label i) const
{
    return speciesNames_[i];
}

inline Foam::label
Foam::plasmaSpecies::speciesID(const Foam::word& name) const
{
    return speciesIDs_[name];
}

//- Electromagnetics model access
inline Foam::electromagneticsModel&
Foam::plasmaSpecies::em() const
{
    return em_;
}

//- Physical properties access
inline const Foam::List<Foam::scalar>&
Foam::plasmaSpecies::speciesChargeNumbers() const
{
    return speciesChargeNumbers_;
}

inline Foam::scalar
Foam::plasmaSpecies::speciesChargeNumber(Foam::label i) const
{
    return speciesChargeNumbers_[i];
}

inline const Foam::PtrList<Foam::dimensionedScalar>&
Foam::plasmaSpecies::speciesCharges() const
{
    return speciesCharges_;
}

inline const Foam::dimensionedScalar&
Foam::plasmaSpecies::speciesCharge(Foam::label i) const
{
    return speciesCharges_[i];
}

inline const Foam::PtrList<Foam::dimensionedScalar>&
Foam::plasmaSpecies::speciesMasses() const
{
    return speciesMasses_;
}

inline const Foam::dimensionedScalar&
Foam::plasmaSpecies::speciesMass(Foam::label i) const
{
    return speciesMasses_[i];
}

inline const Foam::PtrList<Foam::volScalarField>&
Foam::plasmaSpecies::numberDensities() const
{
    return numberDensities_;
}

inline Foam::PtrList<Foam::volScalarField>&
Foam::plasmaSpecies::numberDensities()
{
    return numberDensities_;
}

inline const Foam::volScalarField&
Foam::plasmaSpecies::numberDensity(Foam::label i) const
{
    return numberDensities_[i];
}

inline Foam::volScalarField&
Foam::plasmaSpecies::numberDensity(Foam::label i)
{
    return numberDensities_[i];
}

inline const Foam::List<Foam::scalar>&
Foam::plasmaSpecies::speciesMinNumberDensities() const
{
    return speciesMinNumberDensities_;
}

inline Foam::scalar
Foam::plasmaSpecies::speciesMinNumberDensity(Foam::label i) const
{
    return speciesMinNumberDensities_[i];
}

//- Dictionary Access
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

inline const Foam::HashTable<Foam::dictionary>&
Foam::plasmaSpecies::speciesDicts() const
{
    return speciesDicts_;
}

inline const Foam::dictionary&
Foam::plasmaSpecies::defaultSpeciesDict() const
{
    return defaultSpeciesDict_;
}

//- Background gas access
inline const Foam::word&
Foam::plasmaSpecies::backgroundName() const
{
    return backgroundName_;
}

inline const Foam::dimensionedScalar&
Foam::plasmaSpecies::backgroundDensity() const
{
    return backgroundDensity_;
}

inline const Foam::volScalarField&
Foam::plasmaSpecies::totalNeutralDensity() const
{
    return totalNeutralDensity_;
}

inline const Foam::dictionary&
Foam::plasmaSpecies::backgroundDict() const
{
    return backgroundDict_;
}

//- Group access
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

//- Group Access (Physics/Transport)
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

inline const Foam::labelList&
Foam::plasmaSpecies::constantTemperatureSpeciesIDs() const
{
    return constantTemperatureSpeciesIDs_;
}

inline const Foam::labelList&
Foam::plasmaSpecies::dynamicTemperatureSpeciesIDs() const
{
    return dynamicTemperatureSpeciesIDs_;
}

inline const Foam::labelList&
Foam::plasmaSpecies::followBackgroundTempSpeciesIDs() const
{
    return followBackgroundTempSpeciesIDs_;
}

inline const Foam::labelList&
Foam::plasmaSpecies::solveEnergySpeciesIDs() const
{
    return solveEnergySpeciesIDs_;
}

inline const Foam::labelList&
Foam::plasmaSpecies::fieldTemperatureSpeciesIDs() const
{
    return fieldTemperatureSpeciesIDs_;
}

inline const Foam::labelList&
Foam::plasmaSpecies::reactingSpeciesIDs() const
{
    return reactingSpeciesIDs_;
}

inline const Foam::labelList&
Foam::plasmaSpecies::nonReactingSpeciesIDs() const
{
    return nonReactingSpeciesIDs_;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
