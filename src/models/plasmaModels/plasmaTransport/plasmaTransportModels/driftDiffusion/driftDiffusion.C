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
    const volVectorField& E,
    const surfaceScalarField& phiE
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
        phiE
    ),
    advection_(dict.lookupOrDefault<word>("advection", "implicit")),
    fluxScheme_(dict.lookupOrDefault<word>("fluxScheme", "standard")),
    isExplicit_(advection_ == "explicit")
{
    if (fluxScheme_ != "standard" && fluxScheme_ != "ScharfetterGummel")
    {
        FatalIOErrorInFunction(dict_)
            << "Species '" << species_.speciesNames()[specieIndex_] << "': "
            << "Invalid fluxScheme '" << fluxScheme_ << "'" << nl
            << "Valid options are: (standard ScharfetterGummel)" << nl
            << exit(FatalIOError);
    }

    if (advection_ == "explicit")
    {
        isExplicit_ = true;
    }
    else if (advection_ != "implicit")
    {
        FatalIOErrorInFunction(dict_)
            << "Species '" << species_.speciesNames()[specieIndex_] << "': "
            << "Invalid advection '" << advection_ << "'" << nl
            << "Valid options are: (implicit explicit)" << nl
            << exit(FatalIOError);
    }

    if (fluxScheme_ == "ScharfetterGummel" && isExplicit_)
    {
        WarningInFunction
            << "Species '" << species.speciesNames()[specieIndex] << "': "
            << "advection 'explicit' is ignored because "
            << "Scharfetter-Gummel requires implicit coupling. "
            << "Forcing implicit mode." << endl;
        
        isExplicit_ = false;
        advection_ = "implicit";
    }

    constructModels();
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

tmp<fvScalarMatrix> driftDiffusion::nEqn() const
{
    const volScalarField& n = species_.numberDensity(specieIndex_);
    const scalar Z = species_.speciesChargeNumber(specieIndex_);

    surfaceScalarField phi
    (
        IOobject("phi_tmp", mesh_.time().timeName(), mesh_),
        (
            mobilityModel_->isUniform()
          ? Z * mobilityModel_->muValue() * phiE_
          : Z * fvc::interpolate(mobilityModel_->mu()()) * phiE_
        )
    );

    tmp<fvScalarMatrix> tEqn(fvm::ddt(n));
    fvScalarMatrix& nEqn = tEqn.ref();

    if (diffusivityModel_->isUniform())
    {
        const dimensionedScalar& D = diffusivityModel_->DValue();
        
        if (fluxScheme_ == "ScharfetterGummel")
        {
            nEqn += fvm::ScharfetterGummel(n, phi, D);
        }
        else
        {
            if (isExplicit_)
            {
                nEqn += fvc::div(phi, n);
            }
            else
            {
                nEqn += fvm::div(phi, n);
            }

            nEqn -= fvm::laplacian(D, n);
        }
    }
    else
    {
        tmp<volScalarField> tD = diffusivityModel_->D();

        if (fluxScheme_ == "ScharfetterGummel")
        {
            nEqn += fvm::ScharfetterGummel(n, phi, tD());
        }
        else
        {
            if (isExplicit_)
            {
                nEqn += fvc::div(phi, n);
            }
            else
            {
                nEqn += fvm::div(phi, n);
            }

            nEqn -= fvm::laplacian(tD(), n);
        }
    }

    return tEqn;
} 

void driftDiffusion::updateParticleFlux(surfaceScalarField& flux) const
{
    const volScalarField& n = species_.numberDensity(specieIndex_);
    const scalar Z = species_.speciesChargeNumber(specieIndex_);

    surfaceScalarField phi
    (
        IOobject("phi_tmp", mesh_.time().timeName(), mesh_),
        (
            mobilityModel_->isUniform()
          ? Z * mobilityModel_->muValue() * phiE_
          : Z * fvc::interpolate(mobilityModel_->mu()()) * phiE_
        )
    );

    if (diffusivityModel_->isUniform())
    {
        const dimensionedScalar& D = diffusivityModel_->DValue();
        if (fluxScheme_ == "ScharfetterGummel")
            flux = fvc::ScharfetterGummel(n, phi, D);
        else
            flux = phi*fvc::interpolate(n) - D*fvc::snGrad(n)*mesh_.magSf();
    }
    else
    {
        tmp<volScalarField> tD = diffusivityModel_->D();
        if (fluxScheme_ == "ScharfetterGummel")
            flux = fvc::ScharfetterGummel(n, phi, tD());
        else
            flux = phi*fvc::interpolate(n) - fvc::interpolate(tD())*fvc::snGrad(n)*mesh_.magSf();
    }
}

tmp<volScalarField> driftDiffusion::electricalConductivity() const
{
    const volScalarField& n = species_.numberDensity(specieIndex_);
    const dimensionedScalar qmag = mag(species_.speciesCharge(specieIndex_));

    if (mobilityModel_->isUniform())
        return qmag * mobilityModel_->muValue() * n;
    
    return qmag * mobilityModel_->mu() * n;
}

tmp<volScalarField> driftDiffusion::diffusiveChargeSource() const
{
    const volScalarField& n = species_.numberDensity(specieIndex_);
    const dimensionedScalar q = species_.speciesCharge(specieIndex_);

    if (diffusivityModel_->isUniform())
        return q * fvc::laplacian(diffusivityModel_->DValue(), n);

    return q * fvc::laplacian(diffusivityModel_->D(), n);
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
