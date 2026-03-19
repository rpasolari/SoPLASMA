/*---------------------------------------------------------------------------*\
  File: gasTemperatureEnergyModel.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::gasTemperatureEnergyModel.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "fvc.H"
#include "fvm.H"

#include "gasTemperatureEnergyModel.H"
#include "plasmaEnergy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(gasTemperatureEnergyModel, 0);
addToRunTimeSelectionTable
(
    plasmaEnergyModel,
    gasTemperatureEnergyModel, 
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

gasTemperatureEnergyModel::gasTemperatureEnergyModel
(
    const word& modelName,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex,
    const volVectorField& E
)
:
    plasmaEnergyModel
    (
        modelName, 
        dict, 
        mesh, 
        species, 
        specieIndex, 
        E     
    )
{}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void gasTemperatureEnergyModel::correct()
{
    // No operations here. Temperature is managed by the background gas
}

tmp<fvScalarMatrix> gasTemperatureEnergyModel::eEqn() const
{
    return nullptr;
}

bool gasTemperatureEnergyModel::isUniform() const
{
    // Query the manager to see if T_gas is a scalar or a field
    return mesh_.lookupObject<plasmaEnergy>("plasmaEnergy").isGasUniform();
}

bool gasTemperatureEnergyModel::isConstant() const
{
    // Query the manager to see if T_gas is being solved or is fixed
    return mesh_.lookupObject<plasmaEnergy>("plasmaEnergy").isGasConstant();
}

tmp<volScalarField> gasTemperatureEnergyModel::T() const
{
    // Look up the manager and return the background temperature field
    return mesh_.lookupObject<plasmaEnergy>("plasmaEnergy").Tgas();
}

const dimensionedScalar& gasTemperatureEnergyModel::Tvalue() const
{
    // Look up the manager and return the background temperature value
    return mesh_.lookupObject<plasmaEnergy>("plasmaEnergy").TgasValue();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
