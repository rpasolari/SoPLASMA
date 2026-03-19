/*---------------------------------------------------------------------------*\
  File: plasmaTransport.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaTransport.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "plasmaTransport.H"
#include "plasmaTransportModel.H"
#include "driftDiffusion.H"

#include "plasmaProfiler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(plasmaTransport, 0);

// * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * *  //

void plasmaTransport::constructModels()
{
    // Loop over species and create a transport model for each one
    for (label i = 0; i < species_.nSpecies(); ++i)
    {
        const word& sName = species_.speciesNames()[i];
        const dictionary& sDict = species_.speciesDict(sName);

        if (!sDict.found("transportModel"))
        {
            FatalIOErrorInFunction(sDict)
                << "Species '" << sName
                << "' is missing required entry 'transportModel' in "
                << species_.dictName() << nl
                << exit(FatalIOError);
        }

        word modelName;
        sDict.lookup("transportModel") >> modelName;

        // Look for each species' transportModelCoeffs
        const dictionary& modelDict = sDict.subDict("transportModelCoeffs");

        // Construct the model using the runtime selection system
        transportModels_.set
        (
            i,
            plasmaTransportModel::New
            (
                modelName,
                modelDict, 
                mesh_, 
                species_, 
                i, 
                E_,
                phiE_
            )
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaTransport::plasmaTransport
(
    plasmaSpecies& species,
    const fvMesh& mesh,
    const volVectorField& E,
    const volScalarField& ePotential
)
:
    regIOobject
    (
        IOobject
        (
            "plasmaTransport",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    species_(species),
    E_(E),
    ePotential_(ePotential),
    transportModels_(species.nSpecies()),
    particleFlux_(),
    k_eff_
    (
        IOobject("k_eff", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
    ),
    phiE_
    (
        IOobject("phiE", mesh.time().timeName(), mesh),
        -fvc::snGrad(ePotential) * mesh.magSf()
    )
{
    constructModels();

    particleFlux_.setSize(species.nSpecies());

    for (const label i : species_.mobileSpeciesIDs())
    {
        // Surface Scalar Field
        particleFlux_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "particleFlux_" + species.speciesNames()[i], 
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar
                (
                    "zero",
                    dimensionSet(0, 0, -1, 0, 0, 0, 0), 
                    0.0
                )
            )
        );
    }
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void plasmaTransport::correct(const bool firstIter, const bool finalIter)
{
    const label nSpecies = species_.nSpecies();
    const label eIdx = species_.speciesID("e");
    const label iIdx = species_.speciesID("pIon");
    
    // Electric field flux calculation
    phiE_ = -fvc::snGrad(ePotential_) * mesh_.magSf();

    // Calculate Electric field magnitude
    volScalarField Emag(mag(E_));

    // Protection against division by zero
    volScalarField safeEmag = max
    (
        Emag, 
        dimensionedScalar("minE", Emag.dimensions(), 1.0)
    );

    // Ionization coefficient (alpha) calculation
    dimensionedScalar E_const("E_const", Emag.dimensions(), 2.73e7);
    dimensionedScalar E_pow_const("Ep", pow(Emag.dimensions(), 3), 4.3666e26);
    dimensionedScalar aScale("as", dimensionSet(0, -1, 0, 0, 0, 0, 0), 1.0);

    volScalarField alpha = 
        (1.1944e6 + E_pow_const / pow(safeEmag, 3)) 
      * exp(-E_const / safeEmag)
      * aScale;
    
    alpha.correctBoundaryConditions();

    // 5. Source term assembly
    dimensionedScalar eta("eta", alpha.dimensions(), 340.75);
    volScalarField alphaEff = alpha - eta;

    const driftDiffusion& eModel = refCast<const driftDiffusion>
                                                       (transportModels_[eIdx]);

    const auto& mobE = eModel.mobility();
    volScalarField veMag
(
    IOobject("veMag", mesh_.time().timeName(), mesh_),
    mesh_,
    dimensionedScalar("zeroU", dimensionSet(0, 1, -1, 0, 0, 0, 0), 0.0)
);

    if (mobE.isUniform())
    {
        veMag = mobE.muValue() * Emag;
    }
    else
    {
        veMag = mobE.mu() * Emag;
    }


    k_eff_ = alphaEff * veMag;
    k_eff_.correctBoundaryConditions();

    volScalarField explicitSource = k_eff_ * species_.numberDensity(eIdx);

    // Solve Continuity Equations
    tmp<fvScalarMatrix> tEqne = transportModels_[eIdx].nEqn();
    fvScalarMatrix& nEqne = tEqne.ref();
    nEqne -= explicitSource; 
    
    tmp<fvScalarMatrix> tEqni = transportModels_[iIdx].nEqn();
    fvScalarMatrix& nEqni = tEqni.ref();
    nEqni -= explicitSource; 

    nEqni.solve();
    nEqne.solve();

    volScalarField& ne = species_.numberDensity(eIdx);
    dimensionedScalar neMin("neMin", ne.dimensions(), 1e10);
    ne = max(ne, neMin);
    ne.correctBoundaryConditions();

    volScalarField& ni = species_.numberDensity(iIdx);
    dimensionedScalar niMin("niMin", ni.dimensions(), 1e10);
    ni = max(ni, niMin);
    ni.correctBoundaryConditions();

    for (const label i : species_.mobileSpeciesIDs())
    {
        transportModels_[i].updateParticleFlux(particleFlux_[i]);
    }
}

tmp<volScalarField> plasmaTransport::electricalConductivity() const
{
    tmp<volScalarField> tSigma
    (
        new volScalarField
        (
            IOobject
            (
                "electricalConductivity",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "zero", 
                dimensionSet(-1, -3, 3, 0, 0, 2, 0), 
                0.0
            )
        )
    );

    // Get reference
    volScalarField& sigma = tSigma.ref();

    const label nSpecies = species_.nSpecies();
    for (label i = 0; i < nSpecies; ++i)
    {
        sigma.ref() += transportModels_[i].electricalConductivity();
    }

    return tSigma;
}

tmp<volScalarField> plasmaTransport::diffusiveChargeSource() const
{
    tmp<volScalarField> tRhoDiff
    (
        new volScalarField
        (
            IOobject
            (
                "diffusiveChargeSource",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar
            (
                "zero", 
                dimensionSet(0, -3, 0, 0, 0, 1, 0), 
                0.0
            )
        )
    );

    // Get reference
    volScalarField& rhoDiff = tRhoDiff.ref();

    const label nSpecies = species_.nSpecies();
    for (label i = 0; i < nSpecies; ++i)
    {
        rhoDiff.ref() += transportModels_[i].diffusiveChargeSource();
    }

    return tRhoDiff;
}

bool plasmaTransport::writeData(Ostream& os) const
{
    return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //












