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
    // Construct mobility model
    if (!dict_.found("mobilityModel"))
    {
        FatalIOErrorInFunction(dict_)
            << "driftDiffusion requires a 'mobilityModel' entry\n"
            << "Please specify mobilityModel <modelName>;\n"
            << exit(FatalIOError);
    }

    word mobilityModelName;
    dict_.lookup("mobilityModel") >> mobilityModelName;

    if (!dict_.found("mobilityCoeffs"))
    {
        FatalIOErrorInFunction(dict_)
            << "Missing sub-dictionary 'mobilityCoeffs' for mobilityModel '"
            << mobilityModelName << "'\n"
            << "driftDiffusion requires mobilityCoeffs { ... }\n"
            << exit(FatalIOError);
    }

    const dictionary& mobilityCoeffs = dict_.subDict("mobilityCoeffs");

    mobilityModel_.reset
    (
        plasmaMobilityModel::New
        (
            mobilityModelName,
            mobilityCoeffs,
            mesh_,
            species_,
            specieIndex_
        )
    );

    // Construct diffusivity model
    if (!dict_.found("diffusivityModel"))
    {
        FatalIOErrorInFunction(dict_)
            << "driftDiffusion requires a 'diffusivityModel' entry\n"
            << "Please specify diffusivityModel <modelName>;\n"
            << exit(FatalIOError);
    }

    word diffusivityModelName;
    dict_.lookup("diffusivityModel") >> diffusivityModelName;

    if (!dict_.found("diffusivityCoeffs"))
    {
        FatalIOErrorInFunction(dict_)
            << "Missing sub-dictionary 'diffusivityCoeffs' for "
            << "diffusivityModel '" << diffusivityModelName << "'\n"
            << "driftDiffusion requires diffusivityCoeffs { ... }\n"
            << exit(FatalIOError);
    }

    const dictionary& diffusivityCoeffs = dict_.subDict("diffusivityCoeffs");

    diffusivityModel_.reset
    (
        plasmaDiffusivityModel::New
        (
            diffusivityModelName,
            diffusivityCoeffs,
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
    const volVectorField& E,
    const volScalarField& ePotential
)
:
    plasmaTransportModel
    (
        modelName, 
        dict, 
        mesh, 
        species, 
        specieIndex, 
        E, 
        ePotential
    ),
    driftVelocity_(
        IOobject(
            "driftVelocity_" + species.speciesNames()[specieIndex],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector
        (
            "driftVelocity_",
            dimensionSet(0, 1, -1, 0, 0, 0, 0),
            vector::zero
        )
    ),

    phiDrift_(
        IOobject(
            "phiDrift_" + species.speciesNames()[specieIndex],
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "phiDrift_",
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
            IOobject::NO_WRITE
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

void driftDiffusion::correct()
{
    mobilityModel_->correct(mobility_);
    diffusivityModel_->correct(diffusivity_);

    scalar chargeNumber = species_.speciesChargeNumber(specieIndex_);

    driftVelocity_ = chargeNumber * mobility_ * E_;

    phiDrift_ = chargeNumber * fvc::interpolate(mobility_) 
              * (-fvc::snGrad(ePotential_))
              * mesh_.magSf();
}

tmp<fvScalarMatrix> driftDiffusion::nEqn() const
{
    const volScalarField& n = species_.numberDensity(specieIndex_);

    tmp<fvScalarMatrix> tEqn(fvm::ddt(n));
    fvScalarMatrix& nEqn = tEqn.ref();

    if (fluxScheme_ == "ScharfetterGummel")
    {
        nEqn += fvm::ScharfetterGummel(n, phiDrift_, diffusivity_);
    }
    else
    {
        if (isExplicit_)
            nEqn += fvc::div(phiDrift_, n);
        else
            nEqn += fvm::div(phiDrift_, n);

        nEqn -= fvm::laplacian(diffusivity_, n);
    }

    return tEqn;
}

tmp<surfaceScalarField> driftDiffusion::particleFlux() const
{
    const volScalarField& n = species_.numberDensity(specieIndex_);

    tmp<surfaceScalarField> tphi
    (
        new surfaceScalarField
        (
            IOobject
            (
                "flux_" + species_.speciesNames()[specieIndex_],
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
        )
    );

    surfaceScalarField& phi = tphi.ref();

    if (fluxScheme_== "ScharfetterGummel")
    {
        phi = fvc::ScharfetterGummel(n, phiDrift_, diffusivity_);
    }
    else
    {
        phi = phiDrift_*fvc::interpolate(n)
                    - fvc::interpolate(diffusivity_)
                    * fvc::snGrad(n)
                    * mesh_.magSf();
    }

    return tphi; 
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
