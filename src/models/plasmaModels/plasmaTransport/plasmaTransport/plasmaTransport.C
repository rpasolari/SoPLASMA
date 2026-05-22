/*---------------------------------------------------------------------------*\
  File: plasmaTransport.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaTransport.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "plasmaTransport.H"
#include "plasmaTransportModel.H"
#include "plasmaWallBC.H"


// Remove these headers later
#include "interpolationTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Runtime Type Information * * * * * * * * * * //

defineTypeNameAndDebug(plasmaTransport, 0);

// * * * * * * * * * * * * * * Private Member Functions * * * * * * * * * *  //

void plasmaTransport::constructTransportModels()
{
    Info<< "Constructing plasma transport models" << endl;

    // Loop over species and create its transport model
    for (label i = 0; i < species_.nSpecies(); ++i)
    {
        const word& sName = species_.speciesName(i);
        const dictionary& sDict = species_.speciesDict(i);

        if (!sDict.found("transportModel"))
        {
            FatalIOErrorInFunction(sDict)
                << "Species '" << sName
                << "' is missing required entry 'transportModel' in "
                << species_.dictName() << nl << exit(FatalIOError);
        }

        const word modelName(sDict.get<word>("transportModel"));
        const word coeffsName(modelName + "Coeffs");

        const dictionary emptyDict;
        const dictionary& modelDict = sDict.found(coeffsName) 
                                    ? sDict.subDict(coeffsName)
                                    : emptyDict;

        // Construct the model using the RTS
        transportModels_.set
        (
            i,
            plasmaTransportModel::New
            (
                modelName,
                modelDict,
                mesh_,
                species_,
                i
            )
        );

        Info << "    " << sName << ": transport model '" << modelName 
             << "' successfully constructed." << endl;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaTransport::plasmaTransport
(
    const fvMesh& mesh,
    plasmaSpecies& species
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
    transportModels_(species.nSpecies()),
    convectiveFlux_(species.nSpecies()),
    diffusiveFlux_(species.nSpecies()),
    particleFlux_(species.nSpecies()),
    k_eff_
    (
        IOobject
        (
            "k_eff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
    )
{
    constructTransportModels();

    particleFlux_.setSize(species.nSpecies());
    convectiveFlux_.setSize(species.nSpecies());
    diffusiveFlux_.setSize(species.nSpecies());

    const dimensionedScalar zeroFlux
    (
        "zero",
        dimensionSet(0, 0, -1, 0, 0, 0, 0),
        0.0
    );

    for (const label i : species_.mobileSpeciesIDs())
    {
        const word& sName = species_.speciesName(i);

        particleFlux_.set
        (
            i,
            new surfaceScalarField
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
            )
        );

        convectiveFlux_.set
        (
            i,
            new surfaceScalarField
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
            )
        );

        diffusiveFlux_.set
        (
            i,
            new surfaceScalarField
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
            )
        );
    }
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void plasmaTransport::correctTransportModels()
{
    // Update transport coefficients for all models
    for (label i = 0; i < species_.nSpecies(); ++i)
    {
        transportModels_[i].correct();
    }

    Info << "Transport models corrected." << endl;
}

// This is for the positive streamer case
void plasmaTransport::solve()
{
    // ── 1. Update transport coefficients ──────────────────────────────────────
    correctTransportModels();

    // ── 2. Species and fields ─────────────────────────────────────────────────
    const label eIdx = species_.electronSpeciesID();
    const label iIdx = species_.speciesID("pIon");

    volScalarField& ne = species_.numberDensity(eIdx);
    volScalarField& ni = species_.numberDensity(iIdx);

    // Calculate Electric field magnitude
    const volScalarField& Emag = species_.em().Emag();

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

    // Source term assembly
    dimensionedScalar eta("eta", alpha.dimensions(), 340.75);
    volScalarField alphaEff = alpha - eta;

    List<autoPtr<fvScalarMatrix>> eqns(species_.nSpecies());
    for (label i = 0; i < species_.nSpecies(); ++i)
        eqns[i].reset(transportModels_[i].nEqn().ptr());

    for (const label i : species_.mobileSpeciesIDs())
    {
        transportModels_[i].updateFluxes
        (
            *eqns[i],
            convectiveFlux_[i],
            diffusiveFlux_[i],
            particleFlux_[i]
        );
    }

    const volVectorField driftFlux = fvc::reconstruct(convectiveFlux_[eIdx]);

    volScalarField explicitSource_ = alphaEff * mag(driftFlux);
    
    k_eff_ = explicitSource_ / species_.numberDensity(eIdx);
    k_eff_.correctBoundaryConditions();

// const volScalarField mu_e
// (
//     dimensionedScalar("mu_ref", dimensionSet(-1, 0, 2, 0, 0, 1, 0), 2.398)
//   * pow(safeEmag / dimensionedScalar("E_ref", dimensionSet(1, 1, -3, 0, 0, -1, 0), 1.0), scalar(-0.26))
// );

// const volScalarField veMag = mu_e * safeEmag;
//     k_eff_ = alphaEff * veMag;
//     k_eff_.correctBoundaryConditions();
//     volScalarField explicitSource_ = k_eff_ * species_.numberDensity(eIdx);


    // Solve Continuity Equations
    *eqns[eIdx] -= explicitSource_;
    *eqns[iIdx] -= explicitSource_;

    eqns[iIdx]->solve();
    eqns[eIdx]->solve();

    ne.correctBoundaryConditions();
    ni.correctBoundaryConditions();

    species_.clampNumberDensities();
}

// // Solve for SDBD
// void plasmaTransport::solve()
// {
//     // ── 1. Update transport coefficients ──────────────────────────────────────
//     correctTransportModels();

//     // ── 2. Species and fields ─────────────────────────────────────────────────
//     const label eIdx = species_.electronSpeciesID();
//     const label pIdx = species_.speciesID("pIon");
//     const label nIdx = species_.speciesID("nIon");

//     volScalarField& ne = species_.numberDensity(eIdx);
//     volScalarField& ni = species_.numberDensity(pIdx);
//     volScalarField& nn = species_.numberDensity(nIdx);

//     // ── 3. E/N and rate coefficients ──────────────────────────────────────────
//     const dimensionedScalar N_gas("N_gas", dimless/pow(dimLength,3), 2.4463e25);
//     const scalar N_val = N_gas.value();

//     const volScalarField safeEmag = max
//     (
//         species_.em().Emag(),
//         dimensionedScalar("minE", species_.em().Emag().dimensions(), 1.0)
//     );

//     static interpolationTable<scalar> tableAlpha, tableKatt, tableKei;
//     static bool tablesLoaded = false;
//     if (!tablesLoaded)
//     {
//         tableAlpha = interpolationTable<scalar>
//             (mesh_.time().constant()/"totalIonizationReducedTownsendCoeffs");
//         tableKatt  = interpolationTable<scalar>
//             (mesh_.time().constant()/"totalAttachmentRate");
//         tableKei   = interpolationTable<scalar>
//             (mesh_.time().constant()/"totalIonElectronRecombinationRate");
//         tablesLoaded = true;
//     }

//     volScalarField alpha
//     (
//         IOobject("alpha", mesh_.time().timeName(), mesh_,
//                  IOobject::NO_READ, IOobject::AUTO_WRITE),
//         mesh_,
//         dimensionedScalar("zero", dimensionSet(0,-1,0,0,0,0,0), 0.0)
//     );
//     volScalarField k_att
//     (
//         IOobject("k_att", mesh_.time().timeName(), mesh_,
//                  IOobject::NO_READ, IOobject::AUTO_WRITE),
//         mesh_,
//         dimensionedScalar("zero", dimensionSet(0,3,-1,0,0,0,0), 0.0)
//     );
//     volScalarField k_ei(k_att);
//     const dimensionedScalar k_ii("k_ii", dimensionSet(0,3,-1,0,0,0,0), 1.7e-12);

//     forAll(ne, cellI)
//     {
//         const scalar enKey = safeEmag[cellI] / N_val;

//         auto clamp = [](scalar x, const interpolationTable<scalar>& t) -> scalar
//         {
//             return max(t.first().first(), min(t.last().first(), x));
//         };

//         alpha[cellI]  = tableAlpha(clamp(enKey, tableAlpha)) * N_val;
//         k_att[cellI]  = tableKatt(clamp(enKey, tableKatt));
//         k_ei[cellI]   = tableKei(clamp(enKey, tableKei));
//     }

//     // ── 4. Build equations for all mobile species ─────────────────────────────
//     // Immobile species have no equation — mobileSpeciesIDs() skips them
//     List<autoPtr<fvScalarMatrix>> eqns(species_.nSpecies());
//     for (const label i : species_.mobileSpeciesIDs())
//         eqns[i].reset(transportModels_[i].nEqn().ptr());

//     // ── 5. Fill old-n fluxes — needed by chemistry below ─────────────────────
//     // updateFluxes is virtual on plasmaTransportModel — no driftDiffusion cast
//     // Called BEFORE solve: fvMatrix::flux() uses current (old) n
//     for (const label i : species_.mobileSpeciesIDs())
//     {
//         transportModels_[i].updateFluxes
//         (
//             *eqns[i],
//             convectiveFlux_[i],
//             diffusiveFlux_[i],
//             particleFlux_[i]
//         );
//     }
//     // ── 6. Townsend ionization from old-n electron fluxes ────────────────────
//     // Use plasmaTransport arrays directly — no driftDiffusion cast

//     const volVectorField driftVec    = fvc::reconstruct(convectiveFlux_[eIdx]);
//     const volVectorField particleVec = fvc::reconstruct(particleFlux_[eIdx]);

//     const dimensionedScalar smallFlux
//         ("small", driftVec.dimensions(), 1e-6);
//     const dimensionedScalar zeroFlux
//         ("zero",  driftVec.dimensions(), 0.0);

//     const volVectorField driftDir    = driftVec / (mag(driftVec) + smallFlux);

//     const volScalarField ionizationFlux = min
//     (
//         max(particleVec & driftDir, zeroFlux),
//         mag(driftVec)
//     );

//     const volScalarField S_iz = alpha * ionizationFlux;

//     // ── 7. Add chemistry to equations (reactive species only) ─────────────────
//     *eqns[eIdx] -= S_iz;
//     *eqns[eIdx] += (k_att * N_gas + k_ei * ni) * ne;

//     *eqns[pIdx] -= S_iz;
//     *eqns[pIdx] += (k_ei * ne + k_ii * nn) * ni;

//     *eqns[nIdx] -= k_att * N_gas * ne;
//     *eqns[nIdx] += k_ii * ni * nn;

//     // ── 8. Solve all mobile species ───────────────────────────────────────────
//     for (const label i : species_.mobileSpeciesIDs())
//         eqns[i]->solve();

//     // ── 9. Clamp densities ────────────────────────────────────────────────────
//     species_.clampNumberDensities();

//     for (const label i : species_.mobileSpeciesIDs())
//     {
//         transportModels_[i].updateFluxes
//         (
//             *eqns[i],
//             convectiveFlux_[i],
//             diffusiveFlux_[i],
//             particleFlux_[i]
//         );
//     }
// }

void plasmaTransport::updateSurfaceCharge()
{
    const scalar dt = mesh_.time().deltaTValue();

    volScalarField& sigma = species_.em().surfCharge();

    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& p = mesh_.boundary()[patchi];

        if (p.coupled()) continue;

        scalarField& sigmaPatch = sigma.boundaryFieldRef()[patchi];

        // Always reset from previous time step value
        sigmaPatch = sigma.oldTime().boundaryField()[patchi];

        const scalarField& magSf = mesh_.magSf().boundaryField()[patchi];

        for (const label i : species_.mobileSpeciesIDs())
        {
            const fvPatchField<scalar>& pField =
                species_.numberDensity(i).boundaryField()[patchi];

            // Use plasmaWallBC interface 
            const plasmaWallBC* pBC =
                dynamic_cast<const plasmaWallBC*>(&pField);

            if (!pBC || !pBC->enableSurfaceCharging()) continue;

            if (!particleFlux_.set(i)) continue;

            sigmaPatch += species_.speciesCharge(i).value()
                       * particleFlux_[i].boundaryField()[patchi]
                       * dt / magSf;
        }

        sigma.correctBoundaryConditions();
    }

    Info << "Surface charge updated." << endl;
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












