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
    wallFlux_(species.nSpecies())
{
    constructTransportModels();

    particleFlux_.setSize(species.nSpecies());
    convectiveFlux_.setSize(species.nSpecies());
    diffusiveFlux_.setSize(species.nSpecies());
    wallFlux_.setSize(species.nSpecies());

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

        wallFlux_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "wallFlux_" + sName,
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

void plasmaTransport::solve()
{
    // ── 1. Update transport coefficients ──────────────────────────────────────
    correctTransportModels();

    // ── 2. Species and fields ─────────────────────────────────────────────────
    const label eIdx = species_.electronSpeciesID();
    const label pIdx = species_.speciesID("pIon");
    const label nIdx = species_.speciesID("nIon");

    volScalarField& ne = species_.numberDensity(eIdx);
    volScalarField& ni = species_.numberDensity(pIdx);
    volScalarField& nn = species_.numberDensity(nIdx);

    // ── 3. E/N and rate coefficients ──────────────────────────────────────────
    const dimensionedScalar N_gas("N_gas", dimless/pow(dimLength,3), 2.4463e25);
    const scalar N_val = N_gas.value();

    const volScalarField safeEmag = max
    (
        species_.em().Emag(),
        dimensionedScalar("minE", species_.em().Emag().dimensions(), 1.0)
    );

    static interpolationTable<scalar> tableAlpha, tableKatt, tableKei;
    static bool tablesLoaded = false;
    if (!tablesLoaded)
    {
        tableAlpha = interpolationTable<scalar>
            (mesh_.time().constant()/"totalIonizationReducedTownsendCoeffs");
        tableKatt  = interpolationTable<scalar>
            (mesh_.time().constant()/"totalAttachmentRate");
        tableKei   = interpolationTable<scalar>
            (mesh_.time().constant()/"totalIonElectronRecombinationRate");
        tablesLoaded = true;
    }

    volScalarField alpha
    (
        IOobject("alpha", mesh_.time().timeName(), mesh_,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh_,
        dimensionedScalar("zero", dimensionSet(0,-1,0,0,0,0,0), 0.0)
    );
    volScalarField k_att
    (
        IOobject("k_att", mesh_.time().timeName(), mesh_,
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh_,
        dimensionedScalar("zero", dimensionSet(0,3,-1,0,0,0,0), 0.0)
    );
    volScalarField k_ei(k_att);
    const dimensionedScalar k_ii("k_ii", dimensionSet(0,3,-1,0,0,0,0), 1.7e-12);

    forAll(ne, cellI)
    {
        const scalar enKey = safeEmag[cellI] / N_val;

        auto clamp = [](scalar x, const interpolationTable<scalar>& t) -> scalar
        {
            return max(t.first().first(), min(t.last().first(), x));
        };

        alpha[cellI]  = tableAlpha(clamp(enKey, tableAlpha)) * N_val;
        k_att[cellI]  = tableKatt(clamp(enKey, tableKatt));
        k_ei[cellI]   = tableKei(clamp(enKey, tableKei));
    }

    // ── 4. Build equations for all mobile species ─────────────────────────────
    // Immobile species have no equation — mobileSpeciesIDs() skips them
    List<autoPtr<fvScalarMatrix>> eqns(species_.nSpecies());
    for (const label i : species_.mobileSpeciesIDs())
        eqns[i].reset(transportModels_[i].nEqn().ptr());

    // ── 5. Fill old-n fluxes — needed by chemistry below ─────────────────────
    // updateFluxes is virtual on plasmaTransportModel — no driftDiffusion cast
    // Called BEFORE solve: fvMatrix::flux() uses current (old) n
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
    // ── 6. Townsend ionization from old-n electron fluxes ────────────────────
    // Use plasmaTransport arrays directly — no driftDiffusion cast

    const volVectorField driftVec    = fvc::reconstruct(convectiveFlux_[eIdx]);
    const volVectorField particleVec = fvc::reconstruct(particleFlux_[eIdx]);

    const dimensionedScalar smallFlux
        ("small", driftVec.dimensions(), 1e-6);
    const dimensionedScalar zeroFlux
        ("zero",  driftVec.dimensions(), 0.0);

    const volVectorField driftDir    = driftVec / (mag(driftVec) + smallFlux);

    const volScalarField ionizationFlux = min
    (
        max(particleVec & driftDir, zeroFlux),
        mag(driftVec)
    );

    const volScalarField S_iz = alpha * ionizationFlux;

    // ── 7. Add chemistry to equations (reactive species only) ─────────────────
    *eqns[eIdx] -= S_iz;
    *eqns[eIdx] += (k_att * N_gas + k_ei * ni) * ne;

    *eqns[pIdx] -= S_iz;
    *eqns[pIdx] += (k_ei * ne + k_ii * nn) * ni;

    *eqns[nIdx] -= k_att * N_gas * ne;
    *eqns[nIdx] += k_ii * ni * nn;

    // ── 8. Solve all mobile species ───────────────────────────────────────────
    for (const label i : species_.mobileSpeciesIDs())
        eqns[i]->solve();

    // ── 9. Clamp densities ────────────────────────────────────────────────────
    species_.clampNumberDensities();

    for (const label i : species_.mobileSpeciesIDs())
    species_.numberDensity(i).correctBoundaryConditions();

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

    Info << "Transport finished" << endl;

// ════════════════════════════════════════════════════════════════════════
    // DEBUG DIAGNOSTICS
    // ════════════════════════════════════════════════════════════════════════
    {
        Pout.precision(4);
        Pout << scientific;

        const scalar kB      = 1.380649e-23;
        const scalar pi_v    = constant::mathematical::pi;
        const scalar vth4_e  = 0.25*Foam::sqrt(8.0*kB*11600.0/(pi_v*9.109e-31));
        const scalar vth4_pi = 0.25*Foam::sqrt(8.0*kB*300.0  /(pi_v*4.789e-26));

        const volScalarField& muE_vol  =
            mesh_.lookupObject<volScalarField>("mu_" + species_.speciesName(eIdx));
        const volScalarField& DE_vol   =
            mesh_.lookupObject<volScalarField>("D_"  + species_.speciesName(eIdx));
        const volScalarField& muP_vol  =
            mesh_.lookupObject<volScalarField>("mu_" + species_.speciesName(pIdx));
        const volScalarField& DP_vol   =
            mesh_.lookupObject<volScalarField>("D_"  + species_.speciesName(pIdx));
        const volScalarField&    redE  = species_.em().reducedE();
        const volScalarField&    Emag  = species_.em().Emag();
        const surfaceScalarField& phiE = species_.em().phiE();

        // ════════════════════════════════════════════════════════════════
        // PART 1 — Interior cell at x≈5e-6, y≈2.75e-5
        // ════════════════════════════════════════════════════════════════
        {
            const point tgt(5e-6, 2.75e-5, 5e-4);
            label bestC = -1;
            scalar bestD = GREAT;
            forAll(mesh_.C(), ci)
            {
                scalar d = mag(mesh_.C()[ci] - tgt);
                if (d < bestD) { bestD = d; bestC = ci; }
            }

            if (bestC >= 0 && bestD < 5e-5)
            {
                const label ci = bestC;

                const scalar ne_c  = ne[ci];
                const scalar ni_c  = ni[ci];
                const scalar nn_c  = nn[ci];
                const scalar muE_c = muE_vol[ci];
                const scalar DE_c  = DE_vol[ci];
                const scalar muP_c = muP_vol[ci];
                const scalar DP_c  = DP_vol[ci];
                const scalar EN_c  = redE[ci];
                const scalar Em_c  = Emag[ci];

                // Chemistry at cell
                const scalar al_c   = alpha[ci];
                const scalar ionF_c = ionizationFlux[ci];
                const scalar Siz_c  = S_iz[ci];
                const scalar katt_c = k_att[ci];
                const scalar kei_c  = k_ei[ci];
                const scalar kii_v  = k_ii.value();

                // Net source terms per species
                const scalar Se_net = Siz_c - katt_c*N_val*ne_c - kei_c*ni_c*ne_c;
                const scalar Si_net = Siz_c - kei_c*ne_c*ni_c - kii_v*nn_c*ni_c;
                const scalar Sn_net = katt_c*N_val*ne_c - kii_v*ni_c*nn_c;

                // Matrix diagonal and source (includes ddt + chemistry)
                const scalar diagE  = eqns[eIdx]->D()()[ci];
                const scalar srcE   = eqns[eIdx]->source()[ci];
                const scalar diagPI = eqns[pIdx]->D()()[ci];
                const scalar srcPI  = eqns[pIdx]->source()[ci];
                const scalar diagNI = eqns[nIdx]->D()()[ci];
                const scalar srcNI  = eqns[nIdx]->source()[ci];

                // Reconstructed electron fluxes (from first updateFluxes, used for ionization)
                const vector dVec_c = driftVec[ci];
                const vector pVec_c = particleVec[ci];
                const vector dDir_c = driftDir[ci];

                Pout
                    << nl
                    << "╔══════════════════════════════════════════════════════════╗" << nl
                    << "║  CELL DIAG  x≈5e-6  y≈2.75e-5                          ║" << nl
                    << "║  ci=" << ci << "  dist=" << bestD << "  Cc=" << mesh_.C()[ci] << nl
                    << "╠══ E / TRANSPORT ══════════════════════════════════════════╣" << nl
                    << "║  |E|=" << Em_c << "  E/N=" << EN_c << nl
                    << "║  mu_e=" << muE_c << "  D_e=" << DE_c << nl
                    << "║  mu_pi=" << muP_c << "  D_pi=" << DP_c << nl
                    << "║  |ud_e|=" << muE_c*Em_c << "  |ud_pi|=" << muP_c*Em_c << "  [m/s]" << nl
                    << "╠══ DENSITIES [m^-3] ═══════════════════════════════════════╣" << nl
                    << "║  ne=" << ne_c << "  ni=" << ni_c << "  nn=" << nn_c << nl
                    << "╠══ ELECTRON FLUXES (pre-solve, used for S_iz) ══════════════╣" << nl
                    << "║  driftVec:    " << dVec_c << "  [m^-2 s^-1]" << nl
                    << "║  particleVec: " << pVec_c << "  [m^-2 s^-1]" << nl
                    << "║  driftDir:    " << dDir_c << nl
                    << "║  ionizFlux:   " << ionF_c << "  [m^-2 s^-1]" << nl
                    << "╠══ CHEMISTRY [m^-3 s^-1] ══════════════════════════════════╣" << nl
                    << "║  alpha=" << al_c << "  [m^-1]" << nl
                    << "║  S_iz=" << Siz_c << nl
                    << "║  k_att*N=" << katt_c*N_val << "  k_ei=" << kei_c << "  k_ii=" << kii_v << nl
                    << "║  e:  S_iz=" << Siz_c << "  -att=" << katt_c*N_val*ne_c
                    <<               "  -ei="  << kei_c*ni_c*ne_c << "  net=" << Se_net << nl
                    << "║  pi: S_iz=" << Siz_c << "  -ei="  << kei_c*ne_c*ni_c
                    <<               "  -ii="  << kii_v*nn_c*ni_c << "  net=" << Si_net << nl
                    << "║  nn: +att=" << katt_c*N_val*ne_c
                    <<               "  -ii="  << kii_v*ni_c*nn_c << "  net=" << Sn_net << nl
                    << "╠══ MATRIX (diag / RHS source) ═════════════════════════════╣" << nl
                    << "║  e:  D=" << diagE  << "  src=" << srcE  << nl
                    << "║  pi: D=" << diagPI << "  src=" << srcPI << nl
                    << "║  nn: D=" << diagNI << "  src=" << srcNI << nl
                    << "╚══════════════════════════════════════════════════════════╝" << nl;
            }
        }

        // ════════════════════════════════════════════════════════════════
        // PARTS 2 & 3 — Boundary face diagnostics
        // ════════════════════════════════════════════════════════════════

auto faceDiag = [&]
(
    const point& tgt,
    const word& diagLabel,
    const word& patchName,
    const scalar maxDist
)
{
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& p = mesh_.boundary()[patchi];
        if (p.name() != patchName) continue;
        if (p.size() == 0) continue;

                label bestF = -1;
                scalar bestD = GREAT;
                forAll(p, facei)
                {
                    const scalar d = mag(p.Cf()[facei] - tgt);
                    if (d < bestD) { bestD = d; bestF = facei; }
                }
                if (bestF < 0 || bestD > maxDist) continue;

                const label fi    = bestF;
                const label celli = p.faceCells()[fi];

                const scalar magSf = mesh_.magSf().boundaryField()[patchi][fi];
                const vector nf    = mesh_.Sf().boundaryField()[patchi][fi] / magSf;
                const scalar phiEf = phiE.boundaryField()[patchi][fi];
                const scalar Ef_n  = phiEf / magSf;

                const scalar EN_f  = redE.boundaryField()[patchi][fi];
                const scalar EN_c  = redE.internalField()[celli];
                const scalar muE_f = muE_vol.boundaryField()[patchi][fi];
                const scalar DE_f  = DE_vol.boundaryField()[patchi][fi];
                const scalar muP_f = muP_vol.boundaryField()[patchi][fi];
                const scalar DP_f  = DP_vol.boundaryField()[patchi][fi];

                const scalar ne_c  = ne.internalField()[celli];
                const scalar ni_c  = ni.internalField()[celli];
                const scalar nn_c  = nn.internalField()[celli];
                const scalar nf_e  = ne.boundaryField()[patchi][fi];
                const scalar nf_pi = ni.boundaryField()[patchi][fi];

                const scalar udE   = (-1.0)*muE_f*Ef_n;
                const scalar udPi  = (+1.0)*muP_f*Ef_n;
                const bool eDtW    = (udE  > 0.0);
                const bool piDtW   = (udPi > 0.0);

                const scalar thE   = vth4_e  * ne_c * magSf;
                const scalar thPi  = vth4_pi * ni_c * magSf;

                const scalar convE  = convectiveFlux_[eIdx].boundaryField()[patchi][fi];
                const scalar diffE  = diffusiveFlux_[eIdx].boundaryField()[patchi][fi];
                const scalar totE   = particleFlux_[eIdx].boundaryField()[patchi][fi];
                const scalar convPi = convectiveFlux_[pIdx].boundaryField()[patchi][fi];
                const scalar diffPi = diffusiveFlux_[pIdx].boundaryField()[patchi][fi];
                const scalar totPi  = particleFlux_[pIdx].boundaryField()[patchi][fi];

                const scalar expE  = eDtW  ? -(thE  + mag(udE) *ne_c*magSf) : -thE;
                const scalar expPi = piDtW ? -(thPi + mag(udPi)*ni_c*magSf) : -thPi;

                Pout
                    << nl
                    << "╔══════════════════════════════════════════════════════════╗" << nl
                    << "║  FACE DIAG  " << diagLabel << nl
                    << "║  patch=" << p.name() << "  fi=" << fi
                    << "  cell=" << celli << "  dist=" << bestD << nl
                    << "║  Cf=" << p.Cf()[fi] << "  nf=" << nf << "  |Sf|=" << magSf << nl
                    << "╠══ E ══════════════════════════════════════════════════════╣" << nl
                    << "║  phiE=" << phiEf << "  Ef_n=" << Ef_n << "  [V/m]" << nl
                    << "║  E/N face=" << EN_f << "  E/N cell=" << EN_c << "  [Vm^2]" << nl
                    << "╠══ TRANSPORT ═══════════════════════════════════════════════╣" << nl
                    << "║  mu_e=" << muE_f << "  D_e=" << DE_f << nl
                    << "║  mu_pi=" << muP_f << "  D_pi=" << DP_f << nl
                    << "╠══ DENSITIES [m^-3] ════════════════════════════════════════╣" << nl
                    << "║  ne  cell=" << ne_c  << "  face=" << nf_e  << nl
                    << "║  ni  cell=" << ni_c  << "  face=" << nf_pi << nl
                    << "║  nn  cell=" << nn_c << nl
                    << "╠══ DRIFT ════════════════════════════════════════════════════╣" << nl
                    << "║  ud_e=" << udE  << "  → " << (eDtW  ? "TOWARD" : "AWAY") << nl
                    << "║  ud_pi=" << udPi << "  → " << (piDtW ? "TOWARD" : "AWAY") << nl
                    << "╠══ THERMAL & EXPECTED [#/s] ════════════════════════════════╣" << nl
                    << "║  vth4_e=" << vth4_e << "  vth4_pi=" << vth4_pi << "  [m/s]" << nl
                    << "║  therm_e=" << thE  << "  therm_pi=" << thPi << nl
                    << "║  exp_e="  << expE  << "  (" << (eDtW  ? "th+drift" : "th only") << ")" << nl
                    << "║  exp_pi=" << expPi << "  (" << (piDtW ? "th+drift" : "th only") << ")" << nl
                    << "╠══ MATRIX FLUXES [#/s] ═════════════════════════════════════╣" << nl
                    << "║  e:  conv=" << convE  << "  diff=" << diffE  << "  tot=" << totE  << nl
                    << "║  pi: conv=" << convPi << "  diff=" << diffPi << "  tot=" << totPi << nl
                    << "╚══════════════════════════════════════════════════════════╝" << nl;
            }
        };

faceDiag(point(1.5e-5, 0.0,     5e-4), "DIELECTRIC   x=1.5e-5 y=0",      "gas_to_dielectric", 1e-4);
faceDiag(point(5e-7,   2.75e-5, 5e-4), "HV_ELECTRODE x=5e-7  y=2.75e-5", "HV_electrode",      1e-4);
    }
    // ════════════════════════════════════════════════════════════════════════
    // END DEBUG
    // ════════════════════════════════════════════════════════════════════════

    // ── 10. Update surface charge from post-solve consistent fluxes ───────────
    updateSurfaceCharge();

}

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












