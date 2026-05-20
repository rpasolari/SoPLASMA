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
#include "interpolationTable.H"
#include "IFstream.H"
#include "fvcScharfetterGummel.H"

#include "plasmaProfiler.H"

#include <iomanip>  // Required for setprecision
#include <fstream>

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
                species_.E(),
                species_.phiE()
            )
        );
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaTransport::plasmaTransport
(
    plasmaSpecies& species,
    const fvMesh& mesh
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
    particleFlux_(),
    particleFluxReconstructed_(),
    particleFluxCalculated_(),
    driftFlux_(),
    driftFluxReconstructed_(),
    driftFluxCalculated_(),
    diffusiveFlux_(),
    diffusiveFluxReconstructed_(),
    diffusiveFluxCalculated_(),
    ionizationFlux_
    (
        IOobject("ionizationFlux", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("zero", dimensionSet(0, -2, -1, 0, 0, 0, 0), 0.0)
    ),
    driftOrCombinedFluxFlag_
    (
        IOobject("driftOrCombinedFluxFlag", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar("flag", dimless, 0.0) // Dimensionless for 0/1 toggle
    ),
    k_eff_
    (
        IOobject("k_eff", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimensionSet(0, 0, -1, 0, 0, 0, 0), 0.0)
    ),
    explicitSource_
    (
        IOobject("explicitSource", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimensionSet(0, -3, -1, 0, 0, 0, 0), 0.0)
    ),
    S_iz_
    (
        IOobject("S_iz", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimensionSet(0, -3, -1, 0, 0, 0, 0), 0.0)
    ),
    S_att_
    (
        IOobject("S_att", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimensionSet(0, -3, -1, 0, 0, 0, 0), 0.0)
    ),
    S_ei_
    (
        IOobject("S_ei", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimensionSet(0, -3, -1, 0, 0, 0, 0), 0.0)
    ),
    S_ii_
    (
        IOobject("S_ii", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimensionSet(0, -3, -1, 0, 0, 0, 0), 0.0)
    ),
    S_e_net_
    (
        IOobject("S_e_net", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimensionSet(0, -3, -1, 0, 0, 0, 0), 0.0)
    ),
    S_p_net_
    (
        IOobject("S_p_net", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimensionSet(0, -3, -1, 0, 0, 0, 0), 0.0)
    ),
    S_n_net_
    (
        IOobject("S_n_net", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimensionSet(0, -3, -1, 0, 0, 0, 0), 0.0)
    ),
    dne_dt_
    (
        IOobject("dne_dt", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimensionSet(0, -3, -1, 0, 0, 0, 0), 0.0)
    ),
    dni_dt_
    (
        IOobject("dni_dt", mesh.time().timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh,
        dimensionedScalar(dimensionSet(0, -3, -1, 0, 0, 0, 0), 0.0)
    ),
    transportModels_(species.nSpecies())
{
    // Print the requested information
    Info << "Total number of species: " << species_.nSpecies() << endl;
    Info << "Number of mobile species: " << species_.mobileSpeciesIDs().size() << endl;

    constructModels();

    particleFlux_.setSize(species.nSpecies());
    particleFluxReconstructed_.setSize(species.nSpecies());
    particleFluxCalculated_.setSize(species.nSpecies());
    driftFlux_.setSize(species.nSpecies());
    driftFluxReconstructed_.setSize(species.nSpecies());
    driftFluxCalculated_.setSize(species.nSpecies());
    diffusiveFlux_.setSize(species.nSpecies());
    diffusiveFluxReconstructed_.setSize(species.nSpecies());
    diffusiveFluxCalculated_.setSize(species.nSpecies());

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
                    IOobject::AUTO_WRITE
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

        particleFluxReconstructed_.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "particleFluxReconstructed_" + species.speciesNames()[i], 
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector
                (
                    "zero",
                    dimensionSet(0, -2, -1, 0, 0, 0, 0), 
                    vector(0, 0, 0)
                )
            )
        );

        particleFluxCalculated_.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "particleFluxCalculated_" + species.speciesNames()[i], 
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector
                (
                    "zero",
                    dimensionSet(0, -2, -1, 0, 0, 0, 0), 
                    vector(0, 0, 0)
                )
            )
        );

        driftFlux_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "driftFlux_" + species.speciesNames()[i], 
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
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

        driftFluxReconstructed_.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "driftFluxReconstructed_" + species.speciesNames()[i], 
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector
                (
                    "zero",
                    dimensionSet(0, -2, -1, 0, 0, 0, 0), 
                    vector(0, 0, 0)
                )
            )
        );

        driftFluxCalculated_.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "driftFluxCalculated_" + species.speciesNames()[i], 
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector
                (
                    "zero",
                    dimensionSet(0, -2, -1, 0, 0, 0, 0), 
                    vector(0, 0, 0)
                )
            )
        );

        diffusiveFlux_.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "diffusiveFlux_" + species.speciesNames()[i], 
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
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

        diffusiveFluxReconstructed_.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "diffusiveFluxReconstructed_" + species.speciesNames()[i], 
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector
                (
                    "zero",
                    dimensionSet(0, -2, -1, 0, 0, 0, 0), 
                    vector(0, 0, 0)
                )
            )
        );

        diffusiveFluxCalculated_.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "diffusiveFluxCalculated_" + species.speciesNames()[i], 
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedVector
                (
                    "zero",
                    dimensionSet(0, -2, -1, 0, 0, 0, 0), 
                    vector(0, 0, 0)
                )
            )
        );
    }
}

// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

// // This is for the glow discharge case
// void Foam::plasmaTransport::correct(const bool firstIter, const bool finalIter)
// {
//     const label eIdx = species_.speciesID("e");
//     const label iIdx = species_.speciesID("pIon");
    
//     const volScalarField& Emag = species_.Emag();
//     const dimensionSet fluxDims(0, -2, -1, 0, 0, 0, 0);   // [m^-2 s^-1]
//     const dimensionSet sourceDims(0, -3, -1, 0, 0, 0, 0); // [m^-3 s^-1]
//     const dimensionSet alphaDims(0, -1, 0, 0, 0, 0, 0);   // [m^-1]

//     // 1. Protected E-field for E/N calculation
//     volScalarField safeEmag = max(Emag, dimensionedScalar("minE", Emag.dimensions(), VSMALL));

//     // 2. Physical Constants
//     // dimensionedScalar n_gas("n_gas", dimensionSet(0, -3, 0, 0, 0, 0, 0), 3.5288e22);
//     dimensionedScalar n_gas("n_gas", dimensionSet(0, -3, 0, 0, 0, 0, 0), 9.65729e21);
//     dimensionedScalar unitTd("unitTd", dimensionSet(1, 4, -3, 0, 0, -1, 0), 1e-21);

//     const driftDiffusion& eModel = refCast<const driftDiffusion>(transportModels_[eIdx]);

//     // --- FIX: EN must be dimensionless for the exp() functions ---
//     // We divide E [V/m] by (n [m^-3] * 1e-21 [V m^2]). Both sides are [V/m].
//     // Result is a dimensionless number representing E/N in Townsends.
//     volScalarField EN = safeEmag / (n_gas * unitTd);

//     // 3. Alpha Calculation [m^-1]
//     // Each coefficient is (Value * n_gas) to get [m^-1]
//     dimensionedScalar c1("c1", alphaDims, 1.1e-22 * 9.65729e21);
//     dimensionedScalar c2("c2", alphaDims, 5.5e-21 * 9.65729e21);
//     dimensionedScalar c3("c3", alphaDims, 3.2e-20 * 9.65729e21);
//     dimensionedScalar c4("c4", alphaDims, 1.5e-20 * 9.65729e21);

//     volScalarField alpha = 
//           c1 * exp(-72.0 / (EN + VSMALL))
//         + c2 * exp(-187.0 / (EN + VSMALL))
//         + c3 * exp(-700.0 / (EN + VSMALL))
//         - c4 * exp(-10000.0 / (EN + VSMALL));
    
//     alpha = max(alpha, dimensionedScalar("0", alphaDims, 0.0));

//     // 4. Flux Calculation [m^-2 s^-1]
//     const volScalarField& ne = species_.numberDensity(eIdx);
//     tmp<volScalarField> tmuE = eModel.mobility().mu();
//     tmp<volScalarField> tDE = eModel.diffusivity().D();

//     // driftTerm: [m^-3] * [m^2/Vs] * [V/m] = [m^-2 s^-1]
//     volVectorField driftTerm = - ne * tmuE() * species_.E();
//     // diffTerm: [m^2/s] * [m^-4] = [m^-2 s^-1]
//     volVectorField diffTerm = tDE() * fvc::grad(ne);
    
//     volVectorField phiE_vec = driftTerm - diffTerm;

//     // 5. Clamping & Magnitude
//     dimensionedVector maxV("maxV", fluxDims, vector(1e25, 1e25, 1e25));
//     phiE_vec = min(max(phiE_vec, -maxV), maxV);
//     volScalarField phiE_mag = mag(phiE_vec);

//     // 6. Source and Frequency Update
//     // alpha [m^-1] * phiE_mag [m^-2 s^-1] = explicitSource_ [m^-3 s^-1]
//     explicitSource_ = alpha * phiE_mag;
    
//     // k_eff = S / n_e [s^-1]
//     k_eff_ = explicitSource_ / (ne + dimensionedScalar("sn", ne.dimensions(), 1e-10));

//     // Stability Clamping
//     explicitSource_ = min(explicitSource_, dimensionedScalar("maxS", sourceDims, 1e32));
//     k_eff_ = min(k_eff_, dimensionedScalar("maxK", dimensionSet(0,0,-1,0,0,0,0), 1e13));

//     // 7. Solve
//     {
//         tmp<fvScalarMatrix> tEqne = transportModels_[eIdx].nEqn();
//         tmp<fvScalarMatrix> tEqni = transportModels_[iIdx].nEqn();
        
//         tEqni.ref() -= explicitSource_; 
//         tEqni.ref().solve();
        
//         tEqne.ref() -= explicitSource_; 
//         tEqne.ref().solve();
//     }

//     // 8. Cleanup
//     for (const label id : {eIdx, iIdx})
//     {
//         volScalarField& n_field = species_.numberDensity(id);
//         n_field = max(n_field, dimensionedScalar("nMin", n_field.dimensions(), 1e5));
//         n_field.correctBoundaryConditions();
//     }

//     for (const label i : species_.mobileSpeciesIDs())
//     {
//         transportModels_[i].updateParticleFlux(particleFlux_[i]);
//     }
// }


// // This is for the positive streamer case
// void plasmaTransport::correct(const bool firstIter, const bool finalIter)
// {
//     const label nSpecies = species_.nSpecies();
//     const label eIdx = species_.speciesID("e");
//     const label iIdx = species_.speciesID("pIon");
    
//     // Electric field flux calculation
//     // phiE_ = -fvc::snGrad(ePotential_) * mesh_.magSf();

//     // Calculate Electric field magnitude
//     const volScalarField& Emag = species_.Emag();

//     // Protection against division by zero
//     volScalarField safeEmag = max
//     (
//         Emag, 
//         dimensionedScalar("minE", Emag.dimensions(), 1.0)
//     );

//     // Ionization coefficient (alpha) calculation
//     dimensionedScalar E_const("E_const", Emag.dimensions(), 2.73e7);
//     dimensionedScalar E_pow_const("Ep", pow(Emag.dimensions(), 3), 4.3666e26);
//     dimensionedScalar aScale("as", dimensionSet(0, -1, 0, 0, 0, 0, 0), 1.0);

//     volScalarField alpha = 
//         (1.1944e6 + E_pow_const / pow(safeEmag, 3)) 
//       * exp(-E_const / safeEmag)
//       * aScale;
    
//     alpha.correctBoundaryConditions();

//     // 5. Source term assembly
//     dimensionedScalar eta("eta", alpha.dimensions(), 340.75);
//     volScalarField alphaEff = alpha - eta;

//     const driftDiffusion& eModel = refCast<const driftDiffusion>
//                                                        (transportModels_[eIdx]);

//     const auto& mobE = eModel.mobility();
//     volScalarField veMag
// (
//     IOobject("veMag", mesh_.time().timeName(), mesh_),
//     mesh_,
//     dimensionedScalar("zeroU", dimensionSet(0, 1, -1, 0, 0, 0, 0), 0.0)
// );

//     if (mobE.isUniform())
//     {
//         veMag = mobE.muValue() * Emag;
//     }
//     else
//     {
//         veMag = mobE.mu() * Emag;
//     }


//     k_eff_ = alphaEff * veMag;
//     k_eff_.correctBoundaryConditions();

//     explicitSource_ = k_eff_ * species_.numberDensity(eIdx);

//     // Solve Continuity Equations
//     tmp<fvScalarMatrix> tEqne = transportModels_[eIdx].nEqn();
//     fvScalarMatrix& nEqne = tEqne.ref();
//     nEqne -= explicitSource_; 
    
//     tmp<fvScalarMatrix> tEqni = transportModels_[iIdx].nEqn();
//     fvScalarMatrix& nEqni = tEqni.ref();
//     nEqni -= explicitSource_; 

//     nEqni.solve();
//     nEqne.solve();
// volScalarField& ne = species_.numberDensity(eIdx);
// volScalarField& ni = species_.numberDensity(iIdx);
//     ne.correctBoundaryConditions();
//     ni.correctBoundaryConditions();

//     for (const label i : species_.mobileSpeciesIDs())
//     {
//         transportModels_[i].updateParticleFlux(particleFlux_[i]);
//     }

//     // volScalarField& ne = species_.numberDensity(eIdx);
//     dimensionedScalar neMin("neMin", ne.dimensions(), 1e5);
//     ne = max(ne, neMin);
//     ne.correctBoundaryConditions();

//     // volScalarField& ni = species_.numberDensity(iIdx);
//     dimensionedScalar niMin("niMin", ni.dimensions(), 1e5);
//     ni = max(ni, niMin);
//     ni.correctBoundaryConditions();
// }


// void Foam::plasmaTransport::correct(const bool firstIter, const bool finalIter)
// {
//             // 1. Setup Species Indices and Fields
//         const label eIdx = species_.electronSpeciesID();
//         const label pIdx = species_.speciesID("pIon");
//         const label nIdx = species_.speciesID("nIon");

//     // This is for nonOrthogonal Corrections    
//     if (firstIter)
//     {

//         volScalarField& ne = species_.numberDensity(eIdx);
//         volScalarField& ni = species_.numberDensity(pIdx);
//         volScalarField& nn = species_.numberDensity(nIdx);
        
//         dimensionedScalar nMinSmall("nMinSmall", ne.dimensions(), 1e5);
//         dimensionedScalar nMinLarge("nMinLarge", ne.dimensions(), 1e11);

//         dimensionedScalar N_gas("N_gas", dimless/pow(dimLength, 3), 2.4463e25); 
//         scalar N_val = N_gas.value();
//         volScalarField safeEmag = max(species_.Emag(), dimensionedScalar("minE", species_.Emag().dimensions(), 1.0));

//         // Interpolation Tables
//         static interpolationTable<scalar> tableAlpha;
//         static interpolationTable<scalar> tableKatt;
//         static interpolationTable<scalar> tableKei;
//         // static interpolationTable<scalar> tableKION;

//         static bool tablesLoaded = false;
//         if (!tablesLoaded)
//         {
//             tableAlpha = interpolationTable<scalar>(mesh().time().constant()/"totalIonizationReducedTownsendCoeffs");
//             tableKatt  = interpolationTable<scalar>(mesh().time().constant()/"totalAttachmentRate");
//             tableKei   = interpolationTable<scalar>(mesh().time().constant()/"totalIonElectronRecombinationRate");
//             // tableKION  = interpolationTable<scalar>(mesh().time().constant()/"totalIonizationRate");
//             tablesLoaded = true;
//         }

//         // Coefficients and Intermediate Fields
//         volScalarField alpha(
//             IOobject("alpha", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
//             mesh(),
//             dimensionedScalar("zero", dimensionSet(0, -1, 0, 0, 0, 0, 0), 0.0)
//         );
//         volScalarField k_att(
//             IOobject("k_att", mesh().time().timeName(), mesh(), IOobject::NO_READ, IOobject::AUTO_WRITE),
//             mesh(),
//             dimensionedScalar("zero", dimensionSet(0, 3, -1, 0, 0, 0, 0), 0.0)
//         );
//         volScalarField k_ei(k_att);
//         volScalarField k_ion(k_att);

//         forAll(ne, cellI)
//         {
//             scalar enKey = safeEmag[cellI] / N_val;
//             scalar enKAlpha = max(tableAlpha.first().first(), min(tableAlpha.last().first(), enKey));
//             scalar enKKatt  = max(tableKatt.first().first(),  min(tableKatt.last().first(),  enKey));
//             scalar enKKei   = max(tableKei.first().first(),   min(tableKei.last().first(),   enKey));
//             // scalar enKKION  = max(tableKION.first().first(),  min(tableKION.last().first(),  enKey));

//             alpha[cellI] = tableAlpha(enKAlpha) * N_val;
//             k_att[cellI] = tableKatt(enKKatt);
//             k_ei[cellI]  = tableKei(enKKei);
//             // k_ion[cellI] = tableKION(enKKION);
//         }

//         dimensionedScalar k_ii("k_ii", dimensionSet(0, 3, -1, 0, 0, 0, 0), 1.7e-12);


//         const driftDiffusion& eModel = refCast<const driftDiffusion>(transportModels_[eIdx]);
//         tmp<volScalarField> tMobE = eModel.mobility().mu();

//         // Inner Corrections
//         for (int innerIter = 0; innerIter < 1; ++innerIter)
//         {
//             for (const label i : species_.mobileSpeciesIDs())
//             {
//                 transportModels_[i].updateParticleFlux(particleFlux_[i]);
//                 transportModels_[i].updateDriftFlux(driftFlux_[i]);
//                 transportModels_[i].updateDiffusiveFlux(diffusiveFlux_[i]);
//             }
 

//             const dimensionedScalar smallFlux("small", driftFluxReconstructed_[eIdx].dimensions(), 1e-6);
//             const dimensionedScalar zeroFlux("zero", driftFluxReconstructed_[eIdx].dimensions(), 0.0);

//             driftFluxReconstructed_[eIdx] = fvc::reconstruct(driftFlux_[eIdx]);
//             driftFluxReconstructed_[pIdx] = fvc::reconstruct(driftFlux_[pIdx]);
//             // driftFluxCalculated_[eIdx] = -ne * tMobE() * species_.E();

//             diffusiveFluxReconstructed_[eIdx] = fvc::reconstruct(diffusiveFlux_[eIdx]);
//             diffusiveFluxReconstructed_[pIdx] = fvc::reconstruct(diffusiveFlux_[pIdx]);
//             // diffusiveFluxCalculated_[eIdx] = -eModel.diffusivity().D() * fvc::grad(ne);

//             particleFluxReconstructed_[eIdx] = driftFluxReconstructed_[eIdx] + diffusiveFluxReconstructed_[eIdx];
//             particleFluxReconstructed_[pIdx] = driftFluxReconstructed_[pIdx] + diffusiveFluxReconstructed_[pIdx];
//             // particleFluxCalculated_[eIdx] = driftFluxCalculated_[eIdx] + diffusiveFluxCalculated_[eIdx];


            
//             // // volVectorField totalEFlux = particleFluxCalculated_[eIdx];
//             // volVectorField totalEFlux = particleFluxReconstructed_[eIdx];
//             // // volVectorField driftFlux  = driftFluxCalculated_[eIdx];
//             // volVectorField driftFlux = driftFluxReconstructed_[eIdx];

//             // ionizationFlux_ = min(mag(totalEFlux), mag(driftFlux));

//             // Ionization flux calculation, directional
//             volVectorField driftDir = driftFluxReconstructed_[eIdx] / (mag(driftFluxReconstructed_[eIdx]) + smallFlux);
//             volScalarField fluxAlongDrift = particleFluxReconstructed_[eIdx] & driftDir;
//             ionizationFlux_ = max(min(fluxAlongDrift, mag(driftFluxReconstructed_[eIdx])),zeroFlux);
//             ionizationFlux_ = min(ionizationFlux_,dimensionedScalar("maxFlux", ionizationFlux_.dimensions(), 1e35));
//             ionizationFlux_ = max(ionizationFlux_, zeroFlux);
//             Info << "ionizationFlux_: min=" << gMin(ionizationFlux_) << " max=" << gMax(ionizationFlux_) << endl;

//             // Just a print
//             if (innerIter == 0)
//             {
//                 const volScalarField diffMag = mag(driftFluxReconstructed_[eIdx]) - mag(fluxAlongDrift);
//                 driftOrCombinedFluxFlag_ = pos(diffMag);
//                 label diffusionCorrected = static_cast<label>(
//                     gSum(driftOrCombinedFluxFlag_.primitiveField())
//                 );
//                 label totalCells = returnReduce(mesh_.nCells(), sumOp<label>());
//                 label driftOnly = totalCells - diffusionCorrected;
//                 Info << "Flux selection summary (0=Drift, 1=Combined):" << endl
//                     << "    - Total cells: " << totalCells << endl
//                     << "    - Drift-only (Flag 0): " << driftOnly << endl
//                     << "    - Diffusion-corrected (Flag 1): " << diffusionCorrected << endl;
//             }



//             // Volumetric Source Terms [#/m^3/s]
//             S_iz_  = alpha * ionizationFlux_;            // e + M -> e + e + M+
//             S_att_ = k_att * N_gas * ne;                 // e + M -> M-
//             S_ei_  = k_ei  * ne * ni;                    // e + M+ -> M
//             S_ii_  = k_ii  * ni * nn;                    // M+ + M- -> 2M

//             // Net per-species production [m^-3 s^-1]
//             S_e_net_ = S_iz_ - S_att_ - S_ei_;
//             S_p_net_ = S_iz_ - S_ei_  - S_ii_;
//             S_n_net_ = S_att_ - S_ii_;

//             // Frequencies (used internally for the implicit treatment)
//             volScalarField nu_iz   = S_iz_ / max(ne, nMinLarge);
//             volScalarField nu_att  = k_att * N_gas;
//             volScalarField nu_ei_e = k_ei * ni;
//             volScalarField nu_ei_i = k_ei * ne;
//             volScalarField nu_ii_i = k_ii * nn;
//             volScalarField nu_ii_n = k_ii * ni;

//             k_eff_ = max
//             (
//                 max(max(max(nu_ii_n, nu_ii_i),max(nu_ei_i, nu_ei_e)),
//                 max(nu_att,nu_iz)),
//                 dimensionedScalar("floor", k_eff_.dimensions(), 1e-6)
//             );

//             // 7. Solve Electron Equation
//             tmp<fvScalarMatrix> tEqne = transportModels_[eIdx].nEqn();
//             tEqne.ref() -= S_iz_;
//             tEqne.ref() += nu_att * ne;
//             tEqne.ref() += nu_ei_e * ne;
            

//             // 8. Solve Positive Ion Equation
//             tmp<fvScalarMatrix> tEqni = transportModels_[pIdx].nEqn();
//             tEqni.ref() -= S_iz_;
//             tEqni.ref() += nu_ei_i * ni;
//             tEqni.ref() += nu_ii_i * ni;
            
            
//             // 9. Solve Negative Ion Equation
//             tmp<fvScalarMatrix> tEqnn = transportModels_[nIdx].nEqn();
//             tEqnn.ref() -= nu_att * ne;
//             tEqnn.ref() += nu_ii_n * nn;

//             tEqne.ref().solve();
//             tEqni.ref().solve();
//             tEqnn.ref().solve();

//         // // --- 1. Calculate Individual Frequencies [1/s] ---
//         // volScalarField nu_att  = k_att * N_gas;
//         // volScalarField nu_iz   = k_ion * N_gas;
//         // volScalarField nu_ei_e = k_ei * ni; // electron loss via e-i recomb
//         // volScalarField nu_ei_i = k_ei * ne; // ion loss via e-i recomb
//         // volScalarField nu_ii_i = k_ii * nn; // pos ion loss via i-i recomb
//         // volScalarField nu_ii_n = k_ii * ni; // neg ion loss via i-i recomb

//         // // --- 2. Calculate Explicit Source Rates [#/m^3/s] ---
//         // // volScalarField S_iz  = nu_iz  * ne;
//         // // volScalarField S_att = nu_att * ne;

//         // // --- 3. Update k_eff_ (Maximum of individual frequencies, not sums) ---
//         // // We use nested max() to find the single most restrictive timescale in the chemistry
//         // k_eff_ = max
//         // (
//         //     max(nu_iz, nu_att),
//         //     max(nu_ei_e, max(nu_ei_i, max(nu_ii_i, nu_ii_n)))
//         // );

//         // // Apply the floor for numerical safety
//         // k_eff_ = max
//         // (
//         //     k_eff_,
//         //     dimensionedScalar("floor", k_eff_.dimensions(), 1e-6)
//         // );

//         // // --- 7. Solve Electron Equation ---

//         // tmp<fvScalarMatrix> tEqne = transportModels_[eIdx].nEqn();
        
//         // tEqne.ref() -= nu_iz  * ne;
//         // tEqne.ref() += nu_att * ne;
//         // tEqne.ref() += nu_ei_e * ne;
    

//         // // --- 8. Solve Positive Ion Equation ---

//         // tmp<fvScalarMatrix> tEqni = transportModels_[pIdx].nEqn();
        
//         // tEqni.ref() -= nu_iz  * ne;
//         // tEqni.ref() += nu_ei_i * ni;
//         // tEqni.ref() += nu_ii_i * ni;
        
        
//         // // --- 9. Solve Negative Ion Equation ---
        
//         // tmp<fvScalarMatrix> tEqnn = transportModels_[nIdx].nEqn();
//         // tEqnn.ref() -= nu_att * ne;
//         // tEqnn.ref() += nu_ii_n * nn;
        
//         // tEqne.ref().solve();
//         // tEqni.ref().solve();
//         // tEqnn.ref().solve();


//         Info << "ne: min=" << gMin(ne) << " max=" << gMax(ne) << endl;
//         Info << "ni: min=" << gMin(ni) << " max=" << gMax(ni) << endl;

//         // Clamp and Boundary Conditions (essential for the next inner iteration)
//         ne = max(nMinLarge, ne);
//         ni = max(nMinLarge, ni);
//         nn = max(nMinSmall, nn);
        
//         ne.correctBoundaryConditions();
//         ni.correctBoundaryConditions();
//         nn.correctBoundaryConditions();
//     }


//     // Final Flux Update
//     for (const label i : species_.mobileSpeciesIDs())
//     {
//         transportModels_[i].updateParticleFlux(particleFlux_[i]);
//         transportModels_[i].updateDriftFlux(driftFlux_[i]);
//         transportModels_[i].updateDiffusiveFlux(diffusiveFlux_[i]);
//     }

// S_iz_  = alpha * ionizationFlux_;       // ionizationFlux_ should also be refreshed if you want full consistency
// S_att_ = k_att * N_gas * ne;
// S_ei_  = k_ei  * ne * ni;
// S_ii_  = k_ii  * ni * nn;

// S_e_net_ = S_iz_ - S_att_ - S_ei_;
// S_p_net_ = S_iz_ - S_ei_  - S_ii_;

// // Now budget closes consistently at new state
// dne_dt_ = S_e_net_ - fvc::div(particleFlux_[eIdx]);
// dni_dt_ = S_p_net_ - fvc::div(particleFlux_[pIdx]);


//         // --- FINAL FLUX DEBUG PRINT ---
//         const fvBoundaryMesh& boundary = mesh_.boundary();

//         for (const label i : species_.mobileSpeciesIDs())
//         {
//             const word& sName = species_.speciesNames()[i];
            
//             // Get reference to the density field for this species (adjust 'n_' to your actual field name)
//             const volScalarField& n = species_.numberDensity(i); 

//             const surfaceScalarField& totalFluxField = particleFlux_[i];
//             const surfaceScalarField& driftFluxField = driftFlux_[i];
//             const surfaceScalarField& diffFluxField = diffusiveFlux_[i];

//             forAll(boundary, patchi)
//             {
//                 const fvsPatchScalarField& totalPatch = totalFluxField.boundaryField()[patchi];
//                 const fvsPatchScalarField& driftPatch = driftFluxField.boundaryField()[patchi];
//                 const fvsPatchScalarField& diffPatch  = diffFluxField.boundaryField()[patchi];

//                 const fvPatch& p = boundary[patchi];
                
//                 forAll(p, facei)
//                 {
//                     label celli = p.faceCells()[facei];

//                     if ((mag(p.Cf()[facei].x() - 1.5e-6) < 1e-8 && mag(p.Cf()[facei].y() - 0.0) < 1e-8)
//                         || celli == 1190180)
//                     {
//                         // Access nf (face value) and np (adjacent cell value)
//                         scalar nf = n.boundaryField()[patchi][facei];
//                         scalar np = n.internalField()[celli];

//                         // 1. Set the precision on the Pout object itself
//                         Pout.precision(2);

//                         // 2. Use the Foam::scientific manipulator
//                         Pout << scientific 
//                             << "--- Final Flux Summary [" << sName << "] Face: " << facei << " (Cell: " << celli << ") ---" << endl
//                             << "  Densities | nf: " << nf << " | np: " << np << " [#/m^3]" << endl
//                             << "  Fluxes    | Drift: " << driftPatch[facei] 
//                                         << " | Diff: " << diffPatch[facei] 
//                                         << " | Total: " << totalPatch[facei] << " [#/s]" << endl
//                             << "---------------------------------------" << endl;
//                     }
//                 }
//             }
//         }

// }
// } 

void Foam::plasmaTransport::correct(const bool firstIter, const bool finalIter)
{
    if (!firstIter) return;  // only act on first nonOrthogonal corrector

    // ── Species indices ──────────────────────────────────────────────────────
    const label eIdx = species_.electronSpeciesID();
    const label pIdx = species_.speciesID("pIon");
    const label nIdx = species_.speciesID("nIon");

    volScalarField& ne = species_.numberDensity(eIdx);
    volScalarField& ni = species_.numberDensity(pIdx);
    volScalarField& nn = species_.numberDensity(nIdx);

    dimensionedScalar nMinSmall("nMinSmall", ne.dimensions(), 1e5);
    dimensionedScalar nMinLarge("nMinLarge", ne.dimensions(), 1e11);

    // ── Gas number density and E/N ───────────────────────────────────────────
    dimensionedScalar N_gas("N_gas", dimless/pow(dimLength,3), 2.4463e25);
    scalar N_val = N_gas.value();
    volScalarField safeEmag = max
    (
        species_.Emag(),
        dimensionedScalar("minE", species_.Emag().dimensions(), 1.0)
    );

    // ── Load interpolation tables once ───────────────────────────────────────
    static interpolationTable<scalar> tableAlpha;
    static interpolationTable<scalar> tableKatt;
    static interpolationTable<scalar> tableKei;
    static bool tablesLoaded = false;
    if (!tablesLoaded)
    {
        tableAlpha = interpolationTable<scalar>
            (mesh().time().constant()/"totalIonizationReducedTownsendCoeffs");
        tableKatt  = interpolationTable<scalar>
            (mesh().time().constant()/"totalAttachmentRate");
        tableKei   = interpolationTable<scalar>
            (mesh().time().constant()/"totalIonElectronRecombinationRate");
        tablesLoaded = true;
    }

    // ── Rate coefficients ────────────────────────────────────────────────────
    volScalarField alpha
    (
        IOobject("alpha", mesh().time().timeName(), mesh(),
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh(),
        dimensionedScalar("zero", dimensionSet(0,-1,0,0,0,0,0), 0.0)
    );
    volScalarField k_att
    (
        IOobject("k_att", mesh().time().timeName(), mesh(),
                 IOobject::NO_READ, IOobject::AUTO_WRITE),
        mesh(),
        dimensionedScalar("zero", dimensionSet(0,3,-1,0,0,0,0), 0.0)
    );
    volScalarField k_ei(k_att);

    forAll(ne, cellI)
    {
        scalar enKey   = safeEmag[cellI] / N_val;
        scalar enAlpha = max(tableAlpha.first().first(),
                             min(tableAlpha.last().first(), enKey));
        scalar enKatt  = max(tableKatt.first().first(),
                             min(tableKatt.last().first(),  enKey));
        scalar enKei   = max(tableKei.first().first(),
                             min(tableKei.last().first(),   enKey));

        alpha[cellI]  = tableAlpha(enAlpha) * N_val;
        k_att[cellI]  = tableKatt(enKatt);
        k_ei[cellI]   = tableKei(enKei);
    }

    dimensionedScalar k_ii("k_ii", dimensionSet(0,3,-1,0,0,0,0), 1.7e-12);

    const driftDiffusion& eModel =
        refCast<const driftDiffusion>(transportModels_[eIdx]);

    // ════════════════════════════════════════════════════════════════════════
    // STEP 1 — Build matrices from current n (fills flux cache with old n)
    //          Read fluxes BEFORE solve → consistent with matrix sources
    // ════════════════════════════════════════════════════════════════════════

    // Reset cache so nEqn() rebuilds with current n
    for (const label i : species_.mobileSpeciesIDs())
    {
        if (isA<driftDiffusion>(transportModels_[i]))
            refCast<driftDiffusion>(transportModels_[i]).correct();
    }

    // Build all matrices — this fills the cache
    tmp<fvScalarMatrix> tEqne = transportModels_[eIdx].nEqn();
    tmp<fvScalarMatrix> tEqni = transportModels_[pIdx].nEqn();
    tmp<fvScalarMatrix> tEqnn = transportModels_[nIdx].nEqn();

    // Read fluxes — consistent with the matrices just built (old n)
    for (const label i : species_.mobileSpeciesIDs())
    {
        transportModels_[i].updateParticleFlux(particleFlux_[i]);
        transportModels_[i].updateDriftFlux(driftFlux_[i]);
        transportModels_[i].updateDiffusiveFlux(diffusiveFlux_[i]);
    }

    // Reconstruct for source calculations
    const dimensionedScalar smallFlux
        ("small", driftFluxReconstructed_[eIdx].dimensions(), 1e-6);
    const dimensionedScalar zeroFlux
        ("zero",  driftFluxReconstructed_[eIdx].dimensions(), 0.0);

    driftFluxReconstructed_[eIdx]     = fvc::reconstruct(driftFlux_[eIdx]);
    driftFluxReconstructed_[pIdx]     = fvc::reconstruct(driftFlux_[pIdx]);
    diffusiveFluxReconstructed_[eIdx] = fvc::reconstruct(diffusiveFlux_[eIdx]);
    diffusiveFluxReconstructed_[pIdx] = fvc::reconstruct(diffusiveFlux_[pIdx]);
    particleFluxReconstructed_[eIdx]  =
        driftFluxReconstructed_[eIdx] + diffusiveFluxReconstructed_[eIdx];
    particleFluxReconstructed_[pIdx]  =
        driftFluxReconstructed_[pIdx] + diffusiveFluxReconstructed_[pIdx];

    // Ionization flux (directional, along drift)
    volVectorField driftDir =
        driftFluxReconstructed_[eIdx]
        / (mag(driftFluxReconstructed_[eIdx]) + smallFlux);

    volScalarField fluxAlongDrift =
        particleFluxReconstructed_[eIdx] & driftDir;

    ionizationFlux_ = max
    (
        min(fluxAlongDrift, mag(driftFluxReconstructed_[eIdx])),
        zeroFlux
    );
    ionizationFlux_ = min
    (
        ionizationFlux_,
        dimensionedScalar("maxFlux", ionizationFlux_.dimensions(), 1e35)
    );

    Info << "ionizationFlux_: min=" << gMin(ionizationFlux_)
         << " max=" << gMax(ionizationFlux_) << endl;

    // ── Source terms (computed from old-n fluxes, consistent with matrices) ─
    S_iz_  = alpha * ionizationFlux_;
    S_att_ = k_att * N_gas * ne;
    S_ei_  = k_ei  * ne * ni;
    S_ii_  = k_ii  * ni * nn;

    S_e_net_ = S_iz_ - S_att_ - S_ei_;
    S_p_net_ = S_iz_ - S_ei_  - S_ii_;
    S_n_net_ = S_att_ - S_ii_;

    // Implicit frequencies
    volScalarField nu_iz   = S_iz_  / max(ne, nMinLarge);
    volScalarField nu_att  = k_att  * N_gas;
    volScalarField nu_ei_e = k_ei   * ni;
    volScalarField nu_ei_i = k_ei   * ne;
    volScalarField nu_ii_i = k_ii   * nn;
    volScalarField nu_ii_n = k_ii   * ni;

    k_eff_ = max
    (
        max(max(max(nu_ii_n, nu_ii_i), max(nu_ei_i, nu_ei_e)),
            max(nu_att, nu_iz)),
        dimensionedScalar("floor", k_eff_.dimensions(), 1e-6)
    );

    // ── Add sources to matrices ──────────────────────────────────────────────
    tEqne.ref() -= S_iz_;
    tEqne.ref() += nu_att  * ne;
    tEqne.ref() += nu_ei_e * ne;

    tEqni.ref() -= S_iz_;
    tEqni.ref() += nu_ei_i * ni;
    tEqni.ref() += nu_ii_i * ni;

    tEqnn.ref() -= nu_att * ne;
    tEqnn.ref() += nu_ii_n * nn;

    // ── Solve ────────────────────────────────────────────────────────────────
    tEqne.ref().solve();
    tEqni.ref().solve();
    tEqnn.ref().solve();

    Info << "ne: min=" << gMin(ne) << " max=" << gMax(ne) << endl;
    Info << "ni: min=" << gMin(ni) << " max=" << gMax(ni) << endl;

    // ── Clamp and correct BCs ────────────────────────────────────────────────
    ne = max(nMinLarge, ne);
    ni = max(nMinLarge, ni);
    nn = max(nMinSmall, nn);

    ne.correctBoundaryConditions();
    ni.correctBoundaryConditions();
    nn.correctBoundaryConditions();

    // ════════════════════════════════════════════════════════════════════════
    // STEP 2 — Recompute fluxes from NEW n → use for surface charge etc.
    //          Reset cache so nEqn() uses updated n fields
    // ════════════════════════════════════════════════════════════════════════

    for (const label i : species_.mobileSpeciesIDs())
    {
        if (isA<driftDiffusion>(transportModels_[i]))
            refCast<driftDiffusion>(transportModels_[i]).correct();
    }

    // Rebuild matrices with new n (fills cache)
    for (const label i : species_.mobileSpeciesIDs())
    {
        tmp<fvScalarMatrix> dummy = transportModels_[i].nEqn();
        transportModels_[i].updateParticleFlux(particleFlux_[i]);
        transportModels_[i].updateDriftFlux(driftFlux_[i]);
        transportModels_[i].updateDiffusiveFlux(diffusiveFlux_[i]);
    }

    // ── Budget diagnostics (now consistent: new n, new fluxes) ───────────────
    // Recompute sources with new n for budget closure
    forAll(ne, cellI)
    {
        scalar enKey   = safeEmag[cellI] / N_val;
        scalar enAlpha = max(tableAlpha.first().first(),
                             min(tableAlpha.last().first(), enKey));
        alpha[cellI] = tableAlpha(enAlpha) * N_val;
    }

    // Recompute ionizationFlux with new n fluxes
    driftFluxReconstructed_[eIdx]     = fvc::reconstruct(driftFlux_[eIdx]);
    driftFluxReconstructed_[pIdx]     = fvc::reconstruct(driftFlux_[pIdx]);
    diffusiveFluxReconstructed_[eIdx] = fvc::reconstruct(diffusiveFlux_[eIdx]);
    diffusiveFluxReconstructed_[pIdx] = fvc::reconstruct(diffusiveFlux_[pIdx]);
    particleFluxReconstructed_[eIdx]  =
        driftFluxReconstructed_[eIdx] + diffusiveFluxReconstructed_[eIdx];

    driftDir = driftFluxReconstructed_[eIdx]
               / (mag(driftFluxReconstructed_[eIdx]) + smallFlux);
    fluxAlongDrift = particleFluxReconstructed_[eIdx] & driftDir;
    ionizationFlux_ = max
    (
        min(fluxAlongDrift, mag(driftFluxReconstructed_[eIdx])),
        zeroFlux
    );

    S_iz_  = alpha * ionizationFlux_;
    S_att_ = k_att * N_gas * ne;
    S_ei_  = k_ei  * ne * ni;
    S_ii_  = k_ii  * ni * nn;

    S_e_net_ = S_iz_ - S_att_ - S_ei_;
    S_p_net_ = S_iz_ - S_ei_  - S_ii_;

    dne_dt_ = S_e_net_ - fvc::div(particleFlux_[eIdx]);
    dni_dt_ = S_p_net_ - fvc::div(particleFlux_[pIdx]);

  



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

        const driftDiffusion& pModel =
            refCast<const driftDiffusion>(transportModels_[pIdx]);

        const volScalarField& Emag_vol = species_.Emag();

        const surfaceScalarField& phiE =
            mesh().lookupObject<surfaceScalarField>("phiE");

        // ════════════════════════════════════════════════════════════════
        // PART 1 — Interior cell at x≈5e-6, y≈2.75e-5
        // ════════════════════════════════════════════════════════════════
        {
            const point tgt(5e-6, 2.75e-5, 5e-4);
            label bestC = -1;
            scalar bestD = GREAT;
            forAll(mesh().C(), ci)
            {
                scalar d = mag(mesh().C()[ci] - tgt);
                if (d < bestD) { bestD = d; bestC = ci; }
            }
                // mu/D at cell from model
tmp<volScalarField> tMuE = eModel.mobility().mu();
tmp<volScalarField> tDE  = eModel.diffusivity().D();
tmp<volScalarField> tMuP = pModel.mobility().mu();
tmp<volScalarField> tDP  = pModel.diffusivity().D();
            if (bestC >= 0 && bestD < 5e-5)
            {
                const label ci = bestC;

                const scalar ne_c  = ne[ci];
                const scalar ni_c  = ni[ci];
                const scalar nn_c  = nn[ci];
                const scalar Em_c  = Emag_vol[ci];
                const scalar EN_c  = Em_c / N_val;

    const scalar muE_c = tMuE()[ci];
    const scalar DE_c  = tDE()[ci];
    const scalar muP_c = tMuP()[ci];
    const scalar DP_c  = tDP()[ci];

                // Chemistry
                const scalar al_c   = alpha[ci];
                const scalar ionF_c = ionizationFlux_[ci];
                const scalar Siz_c  = S_iz_[ci];
                const scalar katt_c = k_att[ci];
                const scalar kei_c  = k_ei[ci];
                const scalar kii_v  = k_ii.value();

                const scalar Se_net = Siz_c - katt_c*N_val*ne_c - kei_c*ni_c*ne_c;
                const scalar Si_net = Siz_c - kei_c*ne_c*ni_c - kii_v*nn_c*ni_c;
                const scalar Sn_net = katt_c*N_val*ne_c - kii_v*ni_c*nn_c;

                // Matrix from STEP 1 (old-n matrices, chemistry already added)
                const scalar diagE  = tEqne().D()()[ci];
                const scalar srcE   = tEqne().source()[ci];
                const scalar diagPI = tEqni().D()()[ci];
                const scalar srcPI  = tEqni().source()[ci];
                const scalar diagNI = tEqnn().D()()[ci];
                const scalar srcNI  = tEqnn().source()[ci];

                // Reconstructed electron fluxes (new n, from STEP 2)
                const vector dVec_c = driftFluxReconstructed_[eIdx][ci];
                const vector pVec_c = particleFluxReconstructed_[eIdx][ci];
                const vector dDir_c = driftDir[ci];

                Pout
                    << nl
                    << "╔══════════════════════════════════════════════════════════╗" << nl
                    << "║  CELL DIAG  x≈5e-6  y≈2.75e-5                          ║" << nl
                    << "║  ci=" << ci << "  dist=" << bestD << "  Cc=" << mesh().C()[ci] << nl
                    << "╠══ E / TRANSPORT ══════════════════════════════════════════╣" << nl
                    << "║  |E|=" << Em_c << "  E/N=" << EN_c << nl
                    << "║  mu_e=" << muE_c << "  D_e=" << DE_c << nl
                    << "║  mu_pi=" << muP_c << "  D_pi=" << DP_c << nl
                    << "║  |ud_e|=" << muE_c*Em_c << "  |ud_pi|=" << muP_c*Em_c << "  [m/s]" << nl
                    << "╠══ DENSITIES [m^-3] ═══════════════════════════════════════╣" << nl
                    << "║  ne=" << ne_c << "  ni=" << ni_c << "  nn=" << nn_c << nl
                    << "╠══ ELECTRON FLUXES (new n, STEP 2) ════════════════════════╣" << nl
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
                    << "╠══ MATRIX (diag / RHS src, STEP 1 old-n) ══════════════════╣" << nl
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
            forAll(mesh().boundary(), patchi)
            {
                const fvPatch& p = mesh().boundary()[patchi];
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

                const scalar magSf = mesh().magSf().boundaryField()[patchi][fi];
                const vector nf    = mesh().Sf().boundaryField()[patchi][fi] / magSf;
                const scalar phiEf = phiE.boundaryField()[patchi][fi];
                const scalar Ef_n  = phiEf / magSf;

                // EN: face from Emag boundary (consistent with mu/D), cell from Emag internal
                const scalar EN_f  = Emag_vol.boundaryField()[patchi][fi] / N_val;
                const scalar EN_c  = Emag_vol.internalField()[celli] / N_val;

                // mu/D at face via muPatch/DPatch
                tmp<fvPatchScalarField> tmuE =
                    fvPatchField<scalar>::New("calculated", p, ne.internalField());
                eModel.mobility().muPatch(tmuE.ref(), patchi);
                const scalar muE_f = tmuE()[fi];

                tmp<fvPatchScalarField> tDE =
                    fvPatchField<scalar>::New("calculated", p, ne.internalField());
                eModel.diffusivity().DPatch(tDE.ref(), patchi);
                const scalar DE_f = tDE()[fi];

                tmp<fvPatchScalarField> tmuP =
                    fvPatchField<scalar>::New("calculated", p, ni.internalField());
                pModel.mobility().muPatch(tmuP.ref(), patchi);
                const scalar muP_f = tmuP()[fi];

                tmp<fvPatchScalarField> tDP =
                    fvPatchField<scalar>::New("calculated", p, ni.internalField());
                pModel.diffusivity().DPatch(tDP.ref(), patchi);
                const scalar DP_f = tDP()[fi];

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

                // driftFlux_ = convective (drift), diffusiveFlux_, particleFlux_ = total
                const scalar driftE  = driftFlux_[eIdx].boundaryField()[patchi][fi];
                const scalar diffE   = diffusiveFlux_[eIdx].boundaryField()[patchi][fi];
                const scalar totE    = particleFlux_[eIdx].boundaryField()[patchi][fi];
                const scalar driftPi = driftFlux_[pIdx].boundaryField()[patchi][fi];
                const scalar diffPi  = diffusiveFlux_[pIdx].boundaryField()[patchi][fi];
                const scalar totPi   = particleFlux_[pIdx].boundaryField()[patchi][fi];

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
                    << "╠══ MATRIX FLUXES [#/s] (drift/diff/particle) ═══════════════╣" << nl
                    << "║  e:  drift=" << driftE  << "  diff=" << diffE  << "  tot=" << totE  << nl
                    << "║  pi: drift=" << driftPi << "  diff=" << diffPi << "  tot=" << totPi << nl
                    << "╚══════════════════════════════════════════════════════════╝" << nl;
            }
        };

        faceDiag(point(1.5e-5, 0.0,     5e-4), "DIELECTRIC   x=1.5e-5 y=0",      "gas_to_dielectric", 1e-4);
        faceDiag(point(5e-7,   2.75e-5, 5e-4), "HV_ELECTRODE x=5e-7  y=2.75e-5", "HV_electrode",      1e-4);
    }
    // ════════════════════════════════════════════════════════════════════════
    // END DEBUG
    // ════════════════════════════════════════════════════════════════════════

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












