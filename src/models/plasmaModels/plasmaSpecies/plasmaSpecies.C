/*---------------------------------------------------------------------------*\
  File: plasmaSpecies.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaSpecies.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "plasmaSpecies.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(plasmaSpecies, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaSpecies::plasmaSpecies
(
    const fvMesh& mesh,
    electromagneticsModel& em
)
:
    IOdictionary
    (
        IOobject
        (
            "plasmaSpeciesProperties",
            mesh.time().constant(),
            mesh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    em_(em),
    nSpecies_(0),
    speciesNames_(),
    speciesIDs_(),
    speciesMasses_(),
    speciesCharges_(),
    speciesChargeNumbers_(),
    numberDensities_(),
    speciesMinNumberDensities_(),
    speciesDicts_(),
    defaultSpeciesDict_(),
    backgroundName_("none"),
    backgroundDensity_
    (
        "backgroundDensity", 
        dimensionSet(0, -3, 0, 0, 0, 0, 0),
        0.0
    ),
    totalNeutralDensity_
    (
        IOobject
        (
            "totalNeutralDensity",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "zero", 
            dimensionSet(0, -3, 0, 0, 0, 0, 0), 
            0.0
        )
    ),
    backgroundDict_(),
    electronSpeciesID_(-1),
    ionSpeciesIDs_(),
    positiveIonSpeciesIDs_(),
    negativeIonSpeciesIDs_(),
    neutralSpeciesIDs_(),
    chargedSpeciesIDs_(),
    mobileSpeciesIDs_(),
    immobileSpeciesIDs_(),
    constantTemperatureSpeciesIDs_(),
    dynamicTemperatureSpeciesIDs_(),
    followBackgroundTempSpeciesIDs_(),
    solveEnergySpeciesIDs_(),
    fieldTemperatureSpeciesIDs_(),
    reactingSpeciesIDs_(),
    nonReactingSpeciesIDs_()
{
    if (!found("backgroundGas"))
    {
        FatalIOErrorInFunction(*this)
            << "Missing required dictionary 'backgroundGas' in "
            << objectPath() << nl << exit(FatalIOError);
    }
    
    const dictionary& bgDict = subDict("backgroundGas");
    backgroundName_ = bgDict.get<word>("name");
    backgroundDensity_.value() = bgDict.get<scalar>("numberDensity");
    backgroundDict_ = bgDict;

    totalNeutralDensity_ == backgroundDensity_;

    const dictionary bgEnergyDict = backgroundDict_.subOrEmptyDict("energy");
    bool solveBg = bgEnergyDict.getOrDefault<bool>("solve", false);
    bool bgIsField = solveBg || !bgEnergyDict.found("T");

    constantTemperatureSpeciesIDs_.clear();
    dynamicTemperatureSpeciesIDs_.clear();
    solveEnergySpeciesIDs_.clear();
    fieldTemperatureSpeciesIDs_.clear();

    if (!found("activeSpecies"))
    {
        FatalIOErrorInFunction(*this)
            << "Required entry 'activeSpecies' is missing in dictionary "
            << objectPath() << nl << exit(FatalIOError);
    }

    // Read species list
    lookup("activeSpecies") >> speciesNames_;
    nSpecies_ = speciesNames_.size();

    speciesChargeNumbers_.setSize(nSpecies_);
    speciesCharges_.setSize(nSpecies_);
    speciesMasses_.setSize(nSpecies_);
    numberDensities_.setSize(nSpecies_);
    speciesMinNumberDensities_.setSize(nSpecies_);
    speciesDicts_.resize(nSpecies_);

    // Read species properties dictionary
    if (!found("speciesProperties"))
    {
        FatalIOErrorInFunction(*this)
            << "Missing required dictionary 'speciesProperties' in "
            << objectPath() << nl << exit(FatalIOError);
    }

    const dictionary& propsDict = subDict("speciesProperties");

    // Read defaultProperties if present
    defaultSpeciesDict_ = propsDict.subOrEmptyDict("defaultProperties");

    // Loop over all species
    for (label i = 0; i < nSpecies_; ++i)
    {
        const word& sName = speciesNames_[i];
        speciesIDs_.insert(sName, i);

        // Each species must have a dictionary under speciesProperties
        if (!propsDict.found(sName))
        {
            FatalIOErrorInFunction(*this)
                << "Species '" << sName << "' is listed in 'activeSpecies' but "
                << "has no sub-dictionary in " << objectPath() << nl
                << exit(FatalIOError);
        }

        // Build merged properties (defaults + overrides)
        dictionary mergedDict(defaultSpeciesDict_);
        mergedDict.merge(propsDict.subDict(sName));
        speciesDicts_.insert(sName, mergedDict);

        // Read charge number (required)
        if (!mergedDict.found("charge"))
        {
            FatalIOErrorInFunction(*this)
                << "Species '" << sName << "' is missing required entry "
                << "'charge' in " << objectPath() << nl
                << exit(FatalIOError);
        }

        if (!mergedDict.found("mass"))
        {
            FatalIOErrorInFunction(*this)
                << "Species '" << sName << "' is missing required entry "
                << "'mass' in " << objectPath() << nl
                << exit(FatalIOError);
        }

        speciesChargeNumbers_[i] = readScalar(mergedDict.lookup("charge"));
        scalar massValue = readScalar(mergedDict.lookup("mass"));

        speciesCharges_.set
        (
            i,
            new dimensionedScalar
            (
                "q_" + sName,
                constant::plasma::eCharge.dimensions(),
                speciesChargeNumbers_[i] * constant::plasma::eCharge.value()
            )
        );

        speciesMasses_.set
        (
            i,
            new dimensionedScalar
            (
                "m_" + sName,
                constant::plasma::eMass.dimensions(),
                massValue
            )
        );

        numberDensities_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "n_" + sName,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );

        speciesMinNumberDensities_[i] = 
                       mergedDict.getOrDefault<scalar>("minNumberDensity", 0.0);

        // Groups
        scalar Z = speciesChargeNumbers_[i];
        if (Z == -1 && massValue < 1e-29)
        {
            electronSpeciesID_ = i;
            chargedSpeciesIDs_.append(i);
        }

        // Ions (Charged, but not electrons)
        else if (Z != 0)
        {
            chargedSpeciesIDs_.append(i);
            ionSpeciesIDs_.append(i); 
            if (Z > 0)
                positiveIonSpeciesIDs_.append(i);
            else
                negativeIonSpeciesIDs_.append(i);
        }
        // Active neutrals
        else
        {
            neutralSpeciesIDs_.append(i);
        }
            
        // Mobile vs immobile species
        const word transport =
            mergedDict.getOrDefault<word>("transportModel", "immobile");

        if (transport == "immobile")
        {
            immobileSpeciesIDs_.append(i);
        }
        else
        {
            mobileSpeciesIDs_.append(i);
        }

        // Energy groups
        const word energy =
            mergedDict.getOrDefault<word>("energyModel", "isothermal");

        // Constant vs Dynamic Temperature
        if (energy == "isothermal")
        {
            constantTemperatureSpeciesIDs_.append(i);
        }
        else if (energy == "backgroundGas")
        {
            followBackgroundTempSpeciesIDs_.append(i);
            if (bgIsField)
            {
                fieldTemperatureSpeciesIDs_.append(i);
                dynamicTemperatureSpeciesIDs_.append(i);
            }
            else
            {
                constantTemperatureSpeciesIDs_.append(i);
            }
        }
        else if (energy == "localField")
        {
            fieldTemperatureSpeciesIDs_.append(i);
            dynamicTemperatureSpeciesIDs_.append(i);
        }
        else if (energy == "solveEnergy")
        {
            solveEnergySpeciesIDs_.append(i);
            fieldTemperatureSpeciesIDs_.append(i);
            dynamicTemperatureSpeciesIDs_.append(i);
        }
    }
}  

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void Foam::plasmaSpecies::updateChargeDensity()
{
    em_.chargeDensity() == dimensionedScalar
                                        (em_.chargeDensity().dimensions(), 0.0);

    forAll(chargedSpeciesIDs_, i)
    {
        const label id = chargedSpeciesIDs_[i];
        
        em_.chargeDensity() += numberDensities_[id] * speciesCharges_[id];
    }

    em_.chargeDensity().correctBoundaryConditions();

    Info << "Charge density updated." << endl;
}

void plasmaSpecies::clampNumberDensities()
{
    forAll(numberDensities_, i)
    {
        clampNumberDensity(i);
        numberDensities_[i].correctBoundaryConditions();
    }

    Info << "Species' number densities clamped." << endl;
}

void plasmaSpecies::clampNumberDensity(const label i)
{
    if (speciesMinNumberDensities_[i] > 0.0)
    {
        volScalarField& n = numberDensities_[i];
        const dimensionedScalar nMin
        (
            "nMin",
            n.dimensions(),
            speciesMinNumberDensities_[i]
        );
        n = Foam::max(n, nMin);

        n.correctBoundaryConditions();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
