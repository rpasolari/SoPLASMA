/*---------------------------------------------------------------------------*\
  File: isothermalEnergyModel.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::isothermalEnergyModel.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "fvc.H"
#include "fvm.H"

#include "isothermalEnergyModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(isothermalEnergyModel, 0);
addToRunTimeSelectionTable
(
    plasmaEnergyModel,
    isothermalEnergyModel, 
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

isothermalEnergyModel::isothermalEnergyModel
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
    ),
    T_
    (
        "T_" + species.speciesNames()[specieIndex],
        dimensionSet(0, 0, 0, 1, 0, 0, 0),
        dict.get<scalar>("T")
    )
{}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void isothermalEnergyModel::correct()
{
    // Temperature is constant; no update logic required.
}

tmp<fvScalarMatrix> isothermalEnergyModel::eEqn() const
{
    // No energy equation to solve for an isothermal model.
    return nullptr;
}

bool isothermalEnergyModel::isUniform() const
{
    // In the isothermal model, temperature is always uniform.
    return true;
}

bool isothermalEnergyModel::isConstant() const
{
    // In the isothermal model, temperature is always constant in time.
    return true;
}

tmp<volScalarField> isothermalEnergyModel::T() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            T_.name(),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        T_
    );
}

const dimensionedScalar& isothermalEnergyModel::Tvalue() const
{
    return T_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
