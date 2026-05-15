/*---------------------------------------------------------------------------*\
  File: plasmaPropertyEvaluator.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaPropertyEvaluator.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "plasmaPropertyEvaluator.H"
#include "genericPlasmaPropertyTemplates.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(plasmaPropertyEvaluator, 0);
defineRunTimeSelectionTable(plasmaPropertyEvaluator, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaPropertyEvaluator::plasmaPropertyEvaluator
(
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh,
    const word& propertyName,
    const dimensionSet& dims
)
:
    modelType_(modelType),
    mesh_(mesh),
    propertyName_(propertyName),
    dims_(dims)
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<plasmaPropertyEvaluator> plasmaPropertyEvaluator::New
(
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh,
    const word& propertyName,
    const dimensionSet& dims
)
{
    auto* ctorPtr = dictionaryConstructorTable(modelType);

    if (!ctorPtr)
    {
        FatalIOErrorInFunction(dict)
            << "Unknown plasmaPropertyEvaluator type '" << modelType << "'\n"
            << "Valid evaluators are: "
            << dictionaryConstructorTablePtr_->sortedToc() << nl
            << exit(FatalIOError);
    }

    return autoPtr<plasmaPropertyEvaluator>
    (
        ctorPtr(modelType, dict, mesh, propertyName, dims)
    );
}

// * * * * * * * * * * * * * * * Instantiations * * * * * * * * * * * * * * //

//- Constant Value
typedef ConstantProperty<plasmaPropertyEvaluator> constantEval;
static plasmaPropertyEvaluator::adddictionaryConstructorToTable<constantEval>
    addConstantEvalConstructorToTable_("constant");

//- Tabulated 1D
typedef TabulatedProperty1D<plasmaPropertyEvaluator> tabulatedEval1D;
static plasmaPropertyEvaluator::adddictionaryConstructorToTable<tabulatedEval1D>
    addTabulated1DEvalConstructorToTable_("tabulated1D");

//- Tabulated 2D
typedef TabulatedProperty2D<plasmaPropertyEvaluator> tabulatedEval2D;
static plasmaPropertyEvaluator::adddictionaryConstructorToTable<tabulatedEval2D>
    addTabulated2DEvalConstructorToTable_("tabulated2D");

//- Power Law (mu = A * |X|^b)
typedef PowerLawProperty<plasmaPropertyEvaluator> powerLawEval;
static plasmaPropertyEvaluator::adddictionaryConstructorToTable<powerLawEval>
    addPowerLawEvalConstructorToTable_("powerLaw");

//- Function1 (Coded, Polynomial, etc.)
typedef Function1Property<plasmaPropertyEvaluator> function1Eval;
static plasmaPropertyEvaluator::adddictionaryConstructorToTable<function1Eval>
    addFunction1EvalConstructorToTable_("function1");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
