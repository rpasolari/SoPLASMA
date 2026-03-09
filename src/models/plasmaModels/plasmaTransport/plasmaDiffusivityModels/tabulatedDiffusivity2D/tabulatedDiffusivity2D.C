/*---------------------------------------------------------------------------*\
  File: tabulatedDiffusivity2D.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::tabulatedDiffusivity2D.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "volFields.H"

#include "tabulatedDiffusivity2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(tabulatedDiffusivity2D, 0);
addToRunTimeSelectionTable
(
    plasmaDiffusivityModel, 
    tabulatedDiffusivity2D, 
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tabulatedDiffusivity2D::tabulatedDiffusivity2D
(
    const word& modelName,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex
)
:
    plasmaDiffusivityModel(modelName, dict, mesh, species, specieIndex),
    table_(dict),
    lookupVariables_(dict.get<wordList>("lookupVariables"))
{
    if (dict_.found("file"))
    {
        fileName tableFile = dict_.get<fileName>("file");

        if (!isFile(tableFile))
        {
            FatalIOErrorInFunction(dict)
                << "The specified table file \"" << tableFile
                << "\" does not exist!" << nl
                << "Please ensure the path is correct relative to the "
                << "case directory." << exit(FatalIOError);
        }
    }
    else
    {
        // If the 'file' keyword is missing
        FatalIOErrorInFunction(dict)
            << "Keyword 'file' is missing in diffusivityCoeffs for model '" 
            << modelName << "'." << nl
            << "A tabulated model requires a data file path."
            << exit(FatalIOError);
    }

    // Ensure there are exactly two variables for 2D lookup
    if (lookupVariables_.size() != 2)
    {
        FatalIOErrorInFunction(dict)
            << "tabulatedDiffusivity2D requires exactly 2 lookupVariables, "
            << "but " << lookupVariables_.size() << " were provided."
            << exit(FatalIOError);
    }
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void tabulatedDiffusivity2D::correct(volScalarField& D) const
{
    const volScalarField& fieldX =
        mesh_.lookupObject<volScalarField>(lookupVariables_[0]);
    const volScalarField& fieldY =
        mesh_.lookupObject<volScalarField>(lookupVariables_[1]);

    forAll(D, cellI)
    {
        D[cellI] = table_(fieldX[cellI], fieldY[cellI]); 
    }

    D.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
