/*---------------------------------------------------------------------------*\
  File: plasmaDiffusivityModel.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaDiffusivityModel.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "plasmaDiffusivityModel.H"
#include "genericPlasmaPropertyTemplates.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(plasmaDiffusivityModel, 0);
defineRunTimeSelectionTable(plasmaDiffusivityModel, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaDiffusivityModel::plasmaDiffusivityModel
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

autoPtr<plasmaDiffusivityModel> plasmaDiffusivityModel::New
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
            << "Unknown plasmaDiffusivityModel type '" << modelType << "'\n"
            << "Valid models are: "
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalIOError);
    }

    // Construct and return the model
    return autoPtr<plasmaDiffusivityModel>
    (
        ctorPtr(modelType, dict, mesh, species, specieIndex)
    );
}

// * * * * * * * * * * * * * * * Instantiations * * * * * * * * * * * * * * //

// 1. Constant
typedef ConstantProperty<plasmaDiffusivityModel> 
    constantDiffusivity;
static plasmaDiffusivityModel::
    adddictionaryConstructorToTable<constantDiffusivity>
    addConstantDiffusivityConstructorToTable_("constant");

// 2. Tabulated 1D
typedef TabulatedProperty1D<plasmaDiffusivityModel> 
    tabulatedDiffusivity1D;
static plasmaDiffusivityModel::
    adddictionaryConstructorToTable<tabulatedDiffusivity1D>
    addTabulatedDiffusivity1DConstructorToTable_("tabulated1D");

// 3. Tabulated 2D
typedef TabulatedProperty2D<plasmaDiffusivityModel> 
    tabulatedDiffusivity2D;
static plasmaDiffusivityModel::
    adddictionaryConstructorToTable<tabulatedDiffusivity2D>
    addTabulatedDiffusivity2DConstructorToTable_("tabulated2D");

// 4. Power Law
typedef PowerLawProperty<plasmaDiffusivityModel> 
    powerLawDiffusivity;
static plasmaDiffusivityModel::
    adddictionaryConstructorToTable<powerLawDiffusivity>
    addPowerLawDiffusivityConstructorToTable_("powerLaw");

// 5. Function1
typedef Function1Property<plasmaDiffusivityModel> 
    function1Diffusivity;
static plasmaDiffusivityModel::
    adddictionaryConstructorToTable<function1Diffusivity>
    addFunction1DiffusivityConstructorToTable_("function1");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
