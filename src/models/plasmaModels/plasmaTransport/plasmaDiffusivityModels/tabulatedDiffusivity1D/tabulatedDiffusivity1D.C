/*---------------------------------------------------------------------------*\
  File: tabulatedDiffusivity1D.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::tabulatedDiffusivity1D.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "volFields.H"

#include "tabulatedDiffusivity1D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(tabulatedDiffusivity1D, 0);
addToRunTimeSelectionTable
(
    plasmaDiffusivityModel, 
    tabulatedDiffusivity1D, 
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tabulatedDiffusivity1D::tabulatedDiffusivity1D
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
    lookupVariable_(dict.get<word>("lookupVariable"))
{
    if (dict_.found("file"))
    {
        fileName tableFile = dict_.get<fileName>("file").expand();

        if (!isFile(tableFile))
        {
            FatalIOErrorInFunction(dict_)
                << "The specified table file \"" << tableFile
                << "\" does not exist!" << nl
                << "Please ensure the path is correct relative to the "
                << "case directory." << exit(FatalIOError);
        }
    }
    else
    {
        // If the 'file' keyword is missing
        FatalIOErrorInFunction(dict_)
            << "Keyword 'file' is missing in diffusivityCoeffs for model '" 
            << modelName << "'." << nl
            << "A tabulated model requires a data file path."
            << exit(FatalIOError);
    }

    if (lookupVariable_ != "E")
    {
        FatalIOErrorInFunction(dict)
            << "Unsupported lookupVariable '" << lookupVariable_ << "'." << nl
            << "For tabulatedDiffusivity1D, only 'E' (Electric Field magnitude) "
            << "is supported for now." << exit(FatalIOError);
    }
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void tabulatedDiffusivity1D::correct(volScalarField& D) const
{

    const volVectorField& E =
        mesh_.lookupObject<volVectorField>(lookupVariable_);

    // forAll(D, cellI)
    // {
    //     D[cellI] = table_(mag(E[cellI]));
    // }
    D.primitiveFieldRef() = table_.interpolateValues(mag(E.primitiveField()));

    D.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
