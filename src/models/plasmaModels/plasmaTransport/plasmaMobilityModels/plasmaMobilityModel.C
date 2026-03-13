/*---------------------------------------------------------------------------*\
  File: plasmaMobilityModel.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaMobilityModel.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "plasmaMobilityModel.H"
#include "genericPlasmaPropertyTemplates.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(plasmaMobilityModel, 0);
defineRunTimeSelectionTable(plasmaMobilityModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaMobilityModel::plasmaMobilityModel
(
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex
)
:
    modelType_(modelType),
    mesh_(mesh),
    species_(species),
    specieIndex_(specieIndex),
    dict_(dict)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<plasmaMobilityModel> plasmaMobilityModel::New
(
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex
)
{
    // Lookup constructor using function-call operator
    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        const word& sName = species.speciesNames()[specieIndex];

        FatalIOErrorInFunction(dict)
            << "Species '" << sName << "': "
            << "Unknown plasmaMobilityModel type '" << modelType << "'\n"
            << "Valid models are: "
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalIOError);
    }

    // Construct and return the model
    return autoPtr<plasmaMobilityModel>
    (
        ctorPtr(modelType, dict, mesh, species, specieIndex)
    );
}

// * * * * * * * * * * * * * * * Instantiations * * * * * * * * * * * * * * //

// Constant
typedef ConstantProperty<plasmaMobilityModel> constantMobility;
static plasmaMobilityModel::adddictionaryConstructorToTable<constantMobility>
    addConstantMobilityConstructorToTable_("constant");

// Tabulated 1D
typedef TabulatedProperty1D<plasmaMobilityModel> tabulatedMobility1D;
static plasmaMobilityModel::adddictionaryConstructorToTable<tabulatedMobility1D>
    addTabulatedMobility1DConstructorToTable_("tabulated1D");

// Tabulated 2D
typedef TabulatedProperty2D<plasmaMobilityModel> tabulatedMobility2D;
static plasmaMobilityModel::adddictionaryConstructorToTable<tabulatedMobility2D>
    addTabulatedMobility2DConstructorToTable_("tabulated2D");

// Power Law
typedef PowerLawProperty<plasmaMobilityModel> powerLawMobility;
static plasmaMobilityModel::adddictionaryConstructorToTable<powerLawMobility>
    addPowerLawMobilityConstructorToTable_("powerLaw");

// Function1
typedef Function1Property<plasmaMobilityModel> function1Mobility;
static plasmaMobilityModel::adddictionaryConstructorToTable<function1Mobility>
    addCodedMobilityConstructorToTable_("function1");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
