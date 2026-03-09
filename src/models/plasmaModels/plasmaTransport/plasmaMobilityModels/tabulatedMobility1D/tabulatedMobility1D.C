/*---------------------------------------------------------------------------*\
  File: tabulatedMobility1D.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::tabulatedMobility1D.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "volFields.H"

#include "tabulatedMobility1D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(tabulatedMobility1D, 0);
addToRunTimeSelectionTable
(
    plasmaMobilityModel, 
    tabulatedMobility1D, 
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

tabulatedMobility1D::tabulatedMobility1D
(
    const word& modelName,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex
)
:
    plasmaMobilityModel(modelName, dict, mesh, species, specieIndex),
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
            << "Keyword 'file' is missing in mobilityCoeffs for model '" 
            << modelName << "'." << nl
            << "A tabulated model requires a data file path."
            << exit(FatalIOError);
    }

    if (lookupVariable_ != "E")
    {
        FatalIOErrorInFunction(dict)
            << "Unsupported lookupVariable '" << lookupVariable_ << "'." << nl
            << "For tabulatedMobility1D, only 'E' (Electric Field magnitude) "
            << "is supported for now." << exit(FatalIOError);
    }
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void tabulatedMobility1D::correct(volScalarField& mu) const
{

    const volVectorField& E =
        mesh_.lookupObject<volVectorField>(lookupVariable_);

    // forAll(mu, cellI)
    // {
    //     mu[cellI] = table_(mag(E[cellI]));
    // }
    mu.primitiveFieldRef() = table_.interpolateValues(mag(E.primitiveField()));

    mu.correctBoundaryConditions();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
