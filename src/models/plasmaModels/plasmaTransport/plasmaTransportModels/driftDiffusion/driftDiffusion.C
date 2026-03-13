/*---------------------------------------------------------------------------*\
  File: driftDiffusion.C
  Part of: foamPlasmaToolkit
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
    const word& sName = species_.speciesNames()[specieIndex_];

    // Construct mobility model
    if (!dict_.found("mobility"))
    {
        FatalIOErrorInFunction(dict_)
            << "Species '" << sName << "': missing 'mobility' dict." << nl
            << exit(FatalIOError);
    }

    const dictionary& mobilityDict = dict_.subDict("mobility");

    if (!mobilityDict.found("type"))
    {
        FatalIOErrorInFunction(mobilityDict)
            << "Species '" << sName << "': missing 'type' in mobility." << nl
            << exit(FatalIOError);
    }

    word mobilityType = mobilityDict.get<word>("type");

    mobilityModel_.reset
    (
        plasmaMobilityModel::New
        (
            mobilityType,
            mobilityDict,
            mesh_,
            species_,
            specieIndex_
        )
    );

    // Construct diffusivity model
    if (!dict_.found("diffusivity"))
    {
        FatalIOErrorInFunction(dict_)
            << "Species '" << sName << "': missing 'diffusivity' dict." << nl
            << exit(FatalIOError);
    }

    const dictionary& diffusivityDict = dict_.subDict("diffusivity");

    if (!diffusivityDict.found("type"))
    {
        FatalIOErrorInFunction(diffusivityDict)
            << "Species '" << sName << "': missing 'type' in diffusivity." << nl
            << exit(FatalIOError);
    }

    word diffusivityType = diffusivityDict.get<word>("type");
    diffusivityModel_.reset
    (
        plasmaDiffusivityModel::New
        (
            diffusivityType,
            diffusivityDict,
            mesh_,
            species_,
            specieIndex_
        )
    );

}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

driftDiffusion::driftDiffusion
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
    ),
    driftVelocity_(
        IOobject(
            "driftVelocity_" + species.speciesNames()[specieIndex],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "driftVelocity_",
            dimensionSet(0, 1, -1, 0, 0, 0, 0),
            vector::zero
        )
    ),

    phi_(
        IOobject(
            "phi_" + species.speciesNames()[specieIndex],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "phi_",
            dimensionSet(0, 3, -1, 0, 0, 0, 0),
            0.0
        )
    ),

    mobility_
    (
        IOobject
        (
            "mu_" + species.speciesNames()[specieIndex],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(-1, 0, 2, 0, 0, 1, 0), 0.0)
    ),
    
    diffusivity_
    (
        IOobject
        (
            "D_" + species.speciesNames()[specieIndex],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.0)
    ),
    advectionMode_("implicit"),
    isExplicit_(false)
{
    fluxScheme_ = 
        dict_.lookupOrDefault<word>("driftDiffusionFluxScheme", "standard");

    if (fluxScheme_ != "standard" && fluxScheme_ != "ScharfetterGummel")
    {
        FatalIOErrorInFunction(dict_)
            << "Species '" << species_.speciesNames()[specieIndex_] << "': "
            << "Invalid driftDiffusionFluxScheme '" << fluxScheme_ << "'" << nl
            << "Valid options are: (standard ScharfetterGummel)" << nl
            << exit(FatalIOError);
    }

    advectionMode_ = dict_.lookupOrDefault<word>("advectionMode", "implicit");

    if (advectionMode_ == "explicit")
    {
        isExplicit_ = true;
    }
    else if (advectionMode_ != "implicit")
    {
        FatalIOErrorInFunction(dict_)
            << "Species '" << species_.speciesNames()[specieIndex_] << "': "
            << "Invalid advectionMode '" << advectionMode_ << "'" << nl
            << "Valid options are: (implicit explicit)" << nl
            << exit(FatalIOError);
    }

    if (fluxScheme_ == "ScharfetterGummel" && isExplicit_)
    {
        WarningInFunction
            << "Species '" << species.speciesNames()[specieIndex] << "': "
            << "advectionMode 'explicit' is ignored because "
            << "Scharfetter-Gummel requires implicit coupling. "
            << "Forcing implicit mode." << endl;
        
        isExplicit_ = false;
        advectionMode_ = "implicit";
    }

    constructModels();
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void driftDiffusion::correct(const surfaceScalarField& phiE)
{
    mobilityModel_->correct(mobility_);
    diffusivityModel_->correct(diffusivity_);

    scalar chargeNumber = species_.speciesChargeNumber(specieIndex_);

    driftVelocity_ = chargeNumber * mobility_ * E_;
    driftVelocity_.correctBoundaryConditions();

    phi_ = chargeNumber * fvc::interpolate(mobility_) * phiE;
}

tmp<fvScalarMatrix> driftDiffusion::nEqn() const
{
    const volScalarField& n = species_.numberDensity(specieIndex_);

    tmp<fvScalarMatrix> tEqn(fvm::ddt(n));
    fvScalarMatrix& nEqn = tEqn.ref();

    if (fluxScheme_ == "ScharfetterGummel")
    {
        nEqn += fvm::ScharfetterGummel(n, phi_, diffusivity_);
    }
    else
    {
        if (isExplicit_)
        {
            nEqn += fvc::div(phi_, n);
        }
        else
        {
            nEqn += fvm::div(phi_, n);
        }

        nEqn -= fvm::laplacian(diffusivity_, n);
    }

    return tEqn;
} 

void driftDiffusion::updateParticleFlux(surfaceScalarField& flux) const
{
    const volScalarField& n = species_.numberDensity(specieIndex_);

    if (fluxScheme_ == "ScharfetterGummel")
    {
        flux = fvc::ScharfetterGummel(n, phi_, diffusivity_);
    }
    else
    {
        flux = phi_ * fvc::interpolate(n)
             - fvc::interpolate(diffusivity_) * fvc::snGrad(n) * mesh_.magSf();
    }
}

tmp<surfaceScalarField> driftDiffusion::phi() const
{
    return phi_;
}

const volVectorField& driftDiffusion::driftVelocity() const
{
    return driftVelocity_;
}

const volScalarField* driftDiffusion::diffusivity() const
{
    return &diffusivity_;
}

tmp<volScalarField> driftDiffusion::electricalConductivity() const
{
    return 
        mag(species_.speciesCharge(specieIndex_))
      * mobility_
      * species_.numberDensity(specieIndex_);
}

tmp<volScalarField> driftDiffusion::diffusiveChargeSource() const
{
    return species_.speciesCharge(specieIndex_)
           * fvc::laplacian(diffusivity_, species_.numberDensity(specieIndex_));
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
