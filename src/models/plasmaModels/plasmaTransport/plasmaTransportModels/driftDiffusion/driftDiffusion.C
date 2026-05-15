/*---------------------------------------------------------------------------*\
  File: driftDiffusion.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::driftDiffusion.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "fvc.H"
#include "fvm.H"

#include "driftDiffusion.H"
#include "ScharfetterGummel.H"
#include "fvcScharfetterGummel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(driftDiffusion, 0);
addToRunTimeSelectionTable(plasmaTransportModel, driftDiffusion, dictionary);

// * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * *  //

void driftDiffusion::constructModels()
{
    const word& sName = species_.speciesName(specieIndex_);

    if (!dict_.found("mobility"))
    {
        FatalIOErrorInFunction(dict_)
            << "Species '" << sName
            << "': missing 'mobility' sub-dictionary." << nl
            << exit(FatalIOError);
    }

    const dictionary& mobilityDict = dict_.subDict("mobility");

    if (!mobilityDict.found("type"))
    {
        FatalIOErrorInFunction(mobilityDict)
            << "Species '" << sName
            << "': 'mobility' requires entry 'type'." << nl
            << exit(FatalIOError);
    }

    mobilityModel_ = plasmaMobilityModel::New
    (
        mobilityDict.get<word>("type"),
        mobilityDict,
        mesh_,
        species_,
        specieIndex_
    );

    if (!dict_.found("diffusivity"))
    {
        FatalIOErrorInFunction(dict_)
            << "Species '" << sName
            << "': missing 'diffusivity' sub-dictionary." << nl
            << exit(FatalIOError);
    }

    const dictionary& diffusivityDict = dict_.subDict("diffusivity");

    if (!diffusivityDict.found("type"))
    {
        FatalIOErrorInFunction(diffusivityDict)
            << "Species '" << sName
            << "': 'diffusivity' requires entry 'type'." << nl
            << exit(FatalIOError);
    }

    diffusivityModel_ = plasmaDiffusivityModel::New
    (
        diffusivityDict.get<word>("type"),
        diffusivityDict,
        mesh_,
        species_,
        specieIndex_
    );
}

tmp<surfaceScalarField> driftDiffusion::convectivePhi() const
{
    const scalar Z = species_.speciesChargeNumber(specieIndex_);
    const volScalarField& mu = mobilityModel_->mu();

    return tmp<surfaceScalarField>::New
    (
        IOobject
        (
            "phi_" + species_.speciesName(specieIndex_),
            mesh_.time().timeName(),
            mesh_
        ),
        Z * fvc::interpolate(mu) * species_.em().phiE()
    );
}

void driftDiffusion::initFluxFields()
{
    const word& sName = species_.speciesName(specieIndex_);
    const dimensionedScalar zeroFlux
    (
        "zero",
        dimensionSet(0, 0, -1, 0, 0, 0, 0),
        0.0
    );

    convectiveFlux_ = surfaceScalarField
    (
        IOobject
        (
            "convectiveFlux_" + sName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        zeroFlux
    );

    diffusiveFlux_ = surfaceScalarField
    (
        IOobject
        (
            "diffusiveFlux_" + sName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        zeroFlux
    );

    particleFlux_ = surfaceScalarField
    (
        IOobject
        (
            "particleFlux_" + sName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        zeroFlux
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

driftDiffusion::driftDiffusion
(
    const word& modelName,
    const dictionary& dict,
    const fvMesh& mesh,
    plasmaSpecies& species,
    const label specieIndex
)
:
    plasmaTransportModel(modelName, dict, mesh, species, specieIndex),
    advection_(dict.lookupOrDefault<word>("advection", "implicit")),
    fluxScheme_(dict.lookupOrDefault<word>("fluxScheme", "standard")),
    isExplicit_(advection_ == "explicit"),
    mobilityModel_(nullptr),
    diffusivityModel_(nullptr),
    convectiveFlux_
    (
        IOobject
        (
            "convectiveFlux_" + species.speciesName(specieIndex),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
    ),
    diffusiveFlux_
    (
        IOobject
        (
            "diffusiveFlux_" + species.speciesName(specieIndex),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
    ),
    particleFlux_
    (
        IOobject
        (
            "particleFlux_" + species.speciesName(specieIndex),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
    )
{
    if (fluxScheme_ != "standard" && fluxScheme_ != "ScharfetterGummel")
    {
        FatalIOErrorInFunction(dict_)
            << "Species '" << species_.speciesName(specieIndex_)
            << "': unknown fluxScheme '" << fluxScheme_ << "'." << nl
            << "Valid options: (standard | ScharfetterGummel)" << nl
            << exit(FatalIOError);
    }

    if (advection_ == "explicit")
    {
        isExplicit_ = true;
    }
    else if (advection_ != "implicit")
    {
        FatalIOErrorInFunction(dict_)
            << "Species '" << species_.speciesName(specieIndex_)
            << "': unknown advection '" << advection_ << "'." << nl
            << "Valid options: (implicit | explicit)" << nl
            << exit(FatalIOError);
    }

    if (fluxScheme_ == "ScharfetterGummel" && isExplicit_)
    {
        WarningInFunction
            << "Species '" << species_.speciesName(specieIndex_)
            << "': Scharfetter-Gummel requires implicit advection." << nl
            << "Forcing implicit mode." << endl;

        isExplicit_ = false;
        advection_  = "implicit";
    }

    constructModels();
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void driftDiffusion::correct()
{
    mobilityModel_->correct();
    diffusivityModel_->correct();
}

tmp<fvScalarMatrix> driftDiffusion::nEqn() const
{
    volScalarField& n = species_.numberDensity(specieIndex_);
    const volScalarField& D = diffusivityModel_->D();

    tmp<surfaceScalarField> tPhi = convectivePhi();
    const surfaceScalarField& phi = tPhi();

    tmp<fvScalarMatrix> tEqn = fvm::ddt(n);

    if (fluxScheme_ == "ScharfetterGummel")
    {
        tmp<fvScalarMatrix> tSG = fvm::ScharfetterGummel(n, phi, D);

        particleFlux_    = tSG().flux();
        convectiveFlux_  = phi * fvc::interpolate(n);
        diffusiveFlux_   = particleFlux_ - convectiveFlux_;

        tEqn.ref() += tSG;
    }
    else
    {
        if (isExplicit_)
        {
            convectiveFlux_ = phi * fvc::interpolate(n);
            tEqn.ref() += fvc::div(phi, n);
        }
        else
        {
            tmp<fvScalarMatrix> tDrift = fvm::div(phi, n);
            convectiveFlux_ = tDrift().flux();
            tEqn.ref() += tDrift;
        }

        tmp<fvScalarMatrix> tDiff = fvm::laplacian(D, n);
        diffusiveFlux_ = -tDiff().flux();
        tEqn.ref() -= tDiff;

        particleFlux_ = convectiveFlux_ + diffusiveFlux_;
    }

    return tEqn;
}

void driftDiffusion::updateFluxes
(
    surfaceScalarField& convectiveFlux, 
    surfaceScalarField& diffusiveFlux, 
    surfaceScalarField& particleFlux
) const
{
    convectiveFlux = convectiveFlux_;
    diffusiveFlux = diffusiveFlux_;
    particleFlux = particleFlux_;
}

void driftDiffusion::updateWallFlux(surfaceScalarField& wallFlux) const
{
    wallFlux = dimensionedScalar("zero", particleFlux_.dimensions(), 0.0);

    forAll(wallFlux.boundaryField(), patchi)
    {
        if (!wallFlux.boundaryField()[patchi].coupled())
        {
            wallFlux.boundaryFieldRef()[patchi] =
                particleFlux_.boundaryField()[patchi];
        }
    }
}

tmp<volScalarField> driftDiffusion::electricalConductivity() const
{
    const volScalarField& n = species_.numberDensity(specieIndex_);
    const volScalarField& mu = mobilityModel_->mu();
    const dimensionedScalar qmag = mag(species_.speciesCharge(specieIndex_));

    return qmag * mu * n;
}


tmp<volScalarField> driftDiffusion::diffusiveChargeSource() const
{
    const volScalarField& n = species_.numberDensity(specieIndex_);
    const volScalarField& D = diffusivityModel_->D();
    const dimensionedScalar q = species_.speciesCharge(specieIndex_);

    return q * fvc::laplacian(D, n);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
