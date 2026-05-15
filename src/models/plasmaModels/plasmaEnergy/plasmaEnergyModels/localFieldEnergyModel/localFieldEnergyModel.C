/*---------------------------------------------------------------------------*\
  File: localFieldEnergyModel.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::localFieldEnergyModel.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "fvc.H"
#include "fvm.H"

#include "localFieldEnergyModel.H"
#include "plasmaPropertyEvaluator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(localFieldEnergyModel, 0);
addToRunTimeSelectionTable
(
    plasmaEnergyModel,
    localFieldEnergyModel, 
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

localFieldEnergyModel::localFieldEnergyModel
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
        IOobject
        (
            "T_" + species.speciesNames()[specieIndex],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ, 
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimTemperature, 300.0)
    )
{
    if (!dict.found("temperature"))
    {
        FatalIOErrorInFunction(dict)
            << "Species '" << species.speciesNames()[specieIndex] 
            << "': missing 'temperature' dictionary in energyModelCoeffs."
            << exit(FatalIOError);
    }

    const dictionary& TDict = dict.subDict("temperature");
    const word type = TDict.get<word>("type");

    evaluator_ = plasmaPropertyEvaluator::New
    (
        type,
        TDict,
        mesh,
        species,
        specieIndex,
        "Temperature",
        dimensionSet(0, 0, 0, 1, 0, 0, 0)
    );
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void localFieldEnergyModel::correct()
{
    evaluator_->correct(T_);

    T_.correctBoundaryConditions();
}

tmp<fvScalarMatrix> localFieldEnergyModel::eEqn() const
{
    return nullptr;
}

bool localFieldEnergyModel::isUniform() const
{
    return false;
}

bool localFieldEnergyModel::isConstant() const
{
    return false;
}

tmp<volScalarField> localFieldEnergyModel::T() const
{
    return T_;
}

const dimensionedScalar& localFieldEnergyModel::Tvalue() const
{
    // Explicitly call the base class error to get the species name in the log
    return plasmaEnergyModel::Tvalue();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
