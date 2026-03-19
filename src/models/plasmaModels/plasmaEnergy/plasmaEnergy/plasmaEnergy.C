/*---------------------------------------------------------------------------*\
  File: plasmaEnergy.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaEnergy.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "plasmaEnergy.H"
#include "plasmaEnergyModel.H"

#include "plasmaProfiler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(plasmaEnergy, 0);

// * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * *  //

void plasmaEnergy::constructModels()
{
    // Read the backgroundGas properties
    const dictionary& bgGasDict = species_.backgroundGasDict();
    const dictionary& bgGasEnergyDict = bgGasDict.subDict("energy");

    solveGasEnergy_ = bgGasEnergyDict.get<bool>("solve");
    if (solveGasEnergy_)
    {
        isGasTempField_ = true;
        TgasFieldPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "T_gas", 
                    mesh_.time().timeName(), 
                    mesh_, 
                    IOobject::MUST_READ, 
                    IOobject::AUTO_WRITE
                ),
                mesh_
            )
        );
        TgasValue_.value() = 300.0;
    }
    else
    {
        if (bgGasEnergyDict.found("T"))
        {
            isGasTempField_ = false;
            TgasValue_.value() = bgGasEnergyDict.get<scalar>("T");
            TgasFieldPtr_.reset(nullptr);
        }
        else 
        {
            isGasTempField_ = true;
            TgasFieldPtr_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        "T_gas", 
                        mesh_.time().timeName(), 
                        mesh_, 
                        IOobject::MUST_READ, 
                        IOobject::AUTO_WRITE
                    ),
                    mesh_
                )
            );
            TgasValue_ = dimensionedScalar("Tgas_dummy", dimTemperature, 300.0);
        }
    }

    // Loop over species and create an energy model for each one
    for (label i = 0; i < species_.nSpecies(); ++i)
    {
        const word& sName = species_.speciesNames()[i];
        const dictionary& sDict = species_.speciesDict(sName);

        if (!sDict.found("energyModel"))
        {
            FatalIOErrorInFunction(sDict)
                << "Species '" << sName
                << "' is missing required entry 'energyModel' in "
                << species_.dictName() << nl
                << exit(FatalIOError);
        }

        word modelName;
        sDict.lookup("energyModel") >> modelName;

        // Look for each species' energyModelCoeffs
        const dictionary& modelDict = sDict.subDict("energyModelCoeffs");

        // Construct the model using the runtime selection system
        energyModels_.set
        (
            i,
            plasmaEnergyModel::New
            (
                modelName,
                modelDict, 
                mesh_, 
                species_, 
                i, 
                E_
            )
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaEnergy::plasmaEnergy
(
    plasmaSpecies& species,
    const fvMesh& mesh,
    const volVectorField& E
)
:
    regIOobject
    (
        IOobject
        (
            "plasmaEnergy",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    solveGasEnergy_(false),
    isGasTempField_(false),
    TgasValue_("Tgas", dimTemperature, 300.0),
    mesh_(mesh),
    species_(species),
    E_(E),
    energyModels_(species.nSpecies())
{
    constructModels();
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void plasmaEnergy::correct()
{
    // forAll(energyModels_, i)
    // {
    //     energyModels_[i].correct();
    // }
}

tmp<volScalarField> plasmaEnergy::Tgas() const
{
    if (isGasTempField_)
    {
        return *TgasFieldPtr_;
    }
    else
    {
        return tmp<volScalarField>::New
        (
            IOobject
            (
                "Tgas_tmp", 
                mesh_.time().timeName(), 
                mesh_
            ),
            mesh_,
            TgasValue_
        );
    }
}

const dimensionedScalar& plasmaEnergy::TgasValue() const
{
    if (isGasTempField_)
    {
        FatalErrorInFunction
            << "Requested TgasValue() scalar accessor, but background "
            << "temperature is a spatial field." << nl
            << "Use Tgas() instead." << abort(FatalError);
    }

    return TgasValue_;
}

bool plasmaEnergy::writeData(Ostream& os) const
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

