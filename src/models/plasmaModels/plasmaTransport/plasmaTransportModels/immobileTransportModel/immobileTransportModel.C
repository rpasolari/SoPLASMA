/*---------------------------------------------------------------------------*\
  File: immobileTransportModel.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::immobileTransportModel.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "fvc.H"
#include "fvm.H"

#include "immobileTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(immobileTransportModel, 0);
addToRunTimeSelectionTable
(
    plasmaTransportModel,
    immobileTransportModel, 
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

immobileTransportModel::immobileTransportModel
(
    const word& modelName,
    const dictionary& dict,
    const fvMesh& mesh,
    const plasmaSpecies& species,
    const label specieIndex,
    const volVectorField& E
)
:
    plasmaTransportModel
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

void immobileTransportModel::correct(const surfaceScalarField& phiE)
{
    // Do nothing here    
}

tmp<fvScalarMatrix> immobileTransportModel::nEqn() const
{
    const volScalarField& n = species_.numberDensity(specieIndex_);

    return tmp<fvScalarMatrix>(fvm::ddt(n));
}

void immobileTransportModel::updateParticleFlux(surfaceScalarField& flux) const
{
    // flux *= 0.0;
}

tmp<surfaceScalarField> immobileTransportModel::phi() const
{
    return *zeroSurfaceFieldPtr_;
}

const volVectorField& immobileTransportModel::driftVelocity() const
{
    return *zeroVectorFieldPtr_;
}

tmp<volScalarField> immobileTransportModel::electricalConductivity() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "sigma0", 
            mesh_.time().timeName(), 
            mesh_, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(-1, -3, 3, 0, 0, 2, 0), 0.0)
    );
}

tmp<volScalarField> immobileTransportModel::diffusiveChargeSource() const
{
    return tmp<volScalarField>::New
    (
        IOobject
        (
            "diffSrc0", 
            mesh_.time().timeName(), 
            mesh_, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("0", dimensionSet(0, -3, 0, 0, 0, 1, 0), 0.0)
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
