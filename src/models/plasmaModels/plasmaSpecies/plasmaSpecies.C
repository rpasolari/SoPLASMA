/*---------------------------------------------------------------------------*\
  File: plasmaSpecies.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaSpecies.

  Copyright (C) 2025 Rention Pasolari
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
    const volScalarField& ePotential,
    const volVectorField& E
)
:
    IOdictionary
    (
        IOobject
        (
            "plasmaSpecies",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    ePotential_(ePotential),
    E_(E),
    Emag_
    (
        IOobject
        (
            "Emag", 
            mesh.time().timeName(), 
            mesh, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mag(E_)
    ),
    EN_
    (
        IOobject
        (
            "EN", 
            mesh.time().timeName(), 
            mesh, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(1, 4, -3, 0, 0, -1, 0), 0.0)
    ),
    phiE_
    (
        IOobject
        ("phiE", 
            mesh.time().timeName(), 
            mesh
        ),
        -fvc::snGrad(ePotential) * mesh.magSf()
    ),
    nBackground_("nBackground", dimensionSet(0, -3, 0, 0, 0, 0, 0), 0.0),
    speciesNames_(),
    nSpecies_(0),
    speciesChargeNumber_(),
    speciesCharge_(),
    speciesMass_(),
    numberDensity_(),
    electronSpeciesID_(-1),
    ionSpeciesIDs_(),
    positiveIonSpeciesIDs_(),
    negativeIonSpeciesIDs_(),
    neutralSpeciesIDs_(),
    chargedSpeciesIDs_(),
    mobileSpeciesIDs_(),
    constantTemperatureSpeciesIDs_(),
    solveEnergyEquationSpeciesIDs_()
{
    if (!found("backgroundGas"))
    {
        FatalIOErrorInFunction(*this)
            << "Missing required dictionary 'backgroundGas' in "
            << objectPath() << nl << exit(FatalIOError);
    }
    
    backgroundGasDict_ = subDict("backgroundGas");

    nBackground_ = dimensionedScalar
    (
        "N", 
        dimensionSet(0, -3, 0, 0, 0, 0, 0), 
        backgroundGasDict_.get<scalar>("N")
    );

    const dictionary& bgEnergyDict = backgroundGasDict_.subDict("energy");
    bool solveBg = bgEnergyDict.get<bool>("solve");
    bool bgIsField = false;

    if (solveBg)
    {
        bgIsField = true;
    }
    else
    {
        if (bgEnergyDict.found("T"))
        {
            bgIsField = false;
        }
        else
        {
            bgIsField = true;
        }
    }

    constantTemperatureSpeciesIDs_.clear();
    dynamicTemperatureSpeciesIDs_.clear();
    solveEnergyEquationSpeciesIDs_.clear();
    fieldTemperatureSpeciesIDs_.clear();

    if (!found("activeSpecies"))
    {
        FatalIOErrorInFunction(*this)
            << "Required entry 'activeSpecies' is missing in dictionary "
            << objectPath() << nl
            << exit(FatalIOError);
    }

    // Read species list
    lookup("activeSpecies") >> speciesNames_;
    nSpecies_ = speciesNames_.size();

    speciesID_.clear();
    for (label i = 0; i < nSpecies_; ++i)
    {
        speciesID_.insert(speciesNames_[i], i);
    }

    speciesChargeNumber_.setSize(nSpecies_);
    speciesCharge_.setSize(nSpecies_);
    speciesMass_.setSize(nSpecies_);
    numberDensity_.setSize(nSpecies_);

    // Read species properties dictionary
    if (!found("speciesProperties"))
    {
        FatalIOErrorInFunction(*this)
            << "Missing required dictionary 'speciesProperties' in "
            << objectPath() << nl << exit(FatalIOError);
    }

    const dictionary& propsDict = subDict("speciesProperties");

    // Read defaultProperties if present
    if (propsDict.found("defaultProperties"))
    {
        defaultDict_ = propsDict.subDict("defaultProperties");
    }
    else
    {
        defaultDict_.clear();
    }

    // Prepare storage for merged per-species dicts
    speciesDicts_.clear();
    speciesDicts_.resize(nSpecies_);

    // Loop over all species
    for (label i = 0; i < nSpecies_; ++i)
    {
        const word& sName = speciesNames_[i];

        // Each species must have a dictionary under speciesProperties
        if (!propsDict.found(sName))
        {
            FatalIOErrorInFunction(*this)
                << "Species '" << sName << "' is listed in 'activeSpecies' but "
                << "has no sub-dictionary in " << objectPath() << nl
                << exit(FatalIOError);
        }

        // Build merged properties (defaults + overrides)
        dictionary mergedDict(defaultDict_);
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

        speciesChargeNumber_[i] = readScalar(mergedDict.lookup("charge"));
        scalar massValue = readScalar(mergedDict.lookup("mass"));
        scalar Z = readScalar(mergedDict.lookup("charge"));

        speciesCharge_.set
        (
            i,
            new dimensionedScalar
            (
                "q_" + speciesNames_[i],
                constant::plasma::eCharge.dimensions(),
                speciesChargeNumber_[i] * constant::plasma::eCharge.value()
            )
        );

        speciesMass_.set
        (
            i,
            new dimensionedScalar
            (
                "m_" + speciesNames_[i],
                constant::plasma::eMass.dimensions(),
                massValue
            )
        );

        numberDensity_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "n_" + speciesNames_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );

        // Classify in groups
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
            word transport;
            mergedDict.lookup("transportModel") >> transport;
            if (transport == "immobile")
            {
                immobileSpeciesIDs_.append(i);
            }
            else
            {
                mobileSpeciesIDs_.append(i);
            }

        // Classify in energy groups
            word energy;
            mergedDict.lookup("energyModel") >> energy;

            // Constant vs Dynamic Temperature
            if (energy == "isothermal" || 
                                      (energy == "backgroundGas" && !bgIsField))
            {
                constantTemperatureSpeciesIDs_.append(i);
            }
            else
            {
                dynamicTemperatureSpeciesIDs_.append(i);
            }

            if (energy != "isothermal" && 
                                     !(energy == "backgroundGas" && !bgIsField))
            {
                fieldTemperatureSpeciesIDs_.append(i);
            }

            if (energy == "localEnergy")
            {
                solveEnergyEquationSpeciesIDs_.append(i);
            }

    }
}  

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void plasmaSpecies::updateFields()
{
    Emag_ = mag(E_);
    
    const dimensionedScalar nSmall("s", nBackground_.dimensions(), SMALL);
    EN_ = Emag_ / (nBackground_ + nSmall);

    phiE_ =  -fvc::snGrad(ePotential_) * mesh_.magSf();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
