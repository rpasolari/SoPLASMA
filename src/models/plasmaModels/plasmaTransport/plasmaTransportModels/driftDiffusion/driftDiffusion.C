/*---------------------------------------------------------------------------*\
  File: driftDiffusion.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::driftDiffusion.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "fvc.H"
#include "fvm.H"

#include "driftDiffusion.H"
#include "ScharfetterGummel.H"

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
            << "': missing 'mobility' sub-dictionary in '"
            << dict_.name() << "'." << nl
            << exit(FatalIOError);
    }

    const dictionary& mobilityDict = dict_.subDict("mobility");

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
            << "': missing 'diffusivity' sub-dictionary in '"
            << dict_.name() << "'." << nl
            << exit(FatalIOError);
    }

    const dictionary& diffusivityDict = dict_.subDict("diffusivity");

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
    mobilityModel_(nullptr),
    diffusivityModel_(nullptr),
    fluxScheme_(dict.lookupOrDefault<word>("fluxScheme", "standard")),
    cachedPhi_(),
    convectiveFlux_
    (
        IOobject
        (
            "dd_convectiveFlux_" + species.speciesName(specieIndex),
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
            "dd_diffusiveFlux_" + species.speciesName(specieIndex),
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
            "dd_particleFlux_" + species.speciesName(specieIndex),
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

    constructModels();
}

void driftDiffusion::splitSGFlux
(
    const surfaceScalarField& phi,
    const volScalarField& D,
    const volScalarField& n,
    surfaceScalarField& convFlux,
    surfaceScalarField& diffFlux
) const
{
    auto Bern = [](scalar x) -> scalar
    {
        const scalar ax = mag(x);
        if (ax < 1e-4)
        {
            return 1.0 - 0.5*x + (x*x)/12.0 - pow4(x)/720.0;;
        }
        if (x > 200.0)  return 0.0;
        if (x < -200.0) return -x;
        return x / (Foam::exp(x) - 1.0);
    };

    const surfaceScalarField Df = fvc::interpolate(D);

    const scalarField diffCond = Df.primitiveField()
                                * mesh_.magSf().primitiveField()
                                * mesh_.deltaCoeffs().primitiveField();

    const scalarField& phiI = phi.primitiveField();
    const scalarField& nI = n.primitiveField();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();

    scalarField& convI = convFlux.primitiveFieldRef();
    scalarField& diffI = diffFlux.primitiveFieldRef();

    forAll(phiI, fi)
    {
        convI[fi] = phiI[fi] >= 0
            ? phiI[fi] * nI[own[fi]]
            : phiI[fi] * nI[nei[fi]];

        const scalar absPe = mag(phiI[fi]) / (diffCond[fi] + VSMALL);
        diffI[fi] = diffCond[fi] * Bern(absPe)
                  * (nI[own[fi]] - nI[nei[fi]]);
    }

    forAll(phi.boundaryField(), patchi)
    {
        const scalarField& phiBp  = phi.boundaryField()[patchi];
        const scalarField& nBp    = n.boundaryField()[patchi];
        const labelUList&  fc     = mesh_.boundary()[patchi].faceCells();
        scalarField& convBp       = convFlux.boundaryFieldRef()[patchi];
        scalarField& diffBp       = diffFlux.boundaryFieldRef()[patchi];
        const scalarField& partBp = particleFlux_.boundaryField()[patchi];

        forAll(phiBp, fi)
        {
            convBp[fi] = phiBp[fi] >= 0
                ? phiBp[fi] * nI[fc[fi]]
                : phiBp[fi] * nBp[fi];
        }

        diffBp = partBp - convBp;
    }
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void driftDiffusion::correct()
{
    mobilityModel_->correct();
    diffusivityModel_->correct();

    cachedPhi_.clear();
}

tmp<fvScalarMatrix> driftDiffusion::nEqn() const
{
    volScalarField& n = species_.numberDensity(specieIndex_);
    const volScalarField& D = diffusivityModel_->D();

    cachedPhi_.reset(new surfaceScalarField(convectivePhi()()));
    const surfaceScalarField& phi = *cachedPhi_;

    tmp<fvScalarMatrix> tEqn = fvm::ddt(n);

    if (fluxScheme_ == "ScharfetterGummel")
    {
        Info << "Discretizing transport with SG scheme..." << endl;
        fvScalarMatrix sgMat(fvm::ScharfetterGummel(n, phi, D));

        particleFlux_ = sgMat.flux();
        splitSGFlux(phi, D, n, convectiveFlux_, diffusiveFlux_);

        tEqn.ref() += sgMat;
    }
    else
    {
        Info << "Discretizing transport with standard schemes..." << endl;
        fvScalarMatrix convMat(fvm::div(phi, n));
        convectiveFlux_ = convMat.flux();
        tEqn.ref() += convMat;

        fvScalarMatrix diffMat(fvm::laplacian(D, n));
        diffusiveFlux_ = -diffMat.flux();
        tEqn.ref() -= diffMat;

        particleFlux_ = convectiveFlux_ + diffusiveFlux_;
    }

    return tEqn;
}

void driftDiffusion::updateFluxes
(
    const fvScalarMatrix& nEqnMatrix,
    surfaceScalarField& convectiveFlux, 
    surfaceScalarField& diffusiveFlux, 
    surfaceScalarField& particleFlux
) const
{
    if (!cachedPhi_.valid())
        FatalErrorInFunction
            << "updateFluxes() called before nEqn() — no cached phi." << nl
            << exit(FatalError);

    const volScalarField& n = species_.numberDensity(specieIndex_);
    const volScalarField& D = diffusivityModel_->D();
    const surfaceScalarField& phi = *cachedPhi_;

    particleFlux_ = nEqnMatrix.flux();

    if (fluxScheme_ == "ScharfetterGummel")
    {
        splitSGFlux(phi, D, n, convectiveFlux_, diffusiveFlux_);
    }
    else
    {
        fvScalarMatrix convMat(fvm::div(phi, n));
        convectiveFlux_ = convMat.flux();
        diffusiveFlux_  = particleFlux_ - convectiveFlux_;
    }

    convectiveFlux = convectiveFlux_;
    diffusiveFlux  = diffusiveFlux_;
    particleFlux   = particleFlux_;
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
