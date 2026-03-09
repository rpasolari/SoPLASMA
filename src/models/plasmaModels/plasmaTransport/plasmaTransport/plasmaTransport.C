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

        // Construct the model using the runtime selection system
        transportModels_.set
        (
            i,
            plasmaTransportModel::New
            (
                modelName,
                sDict, 
                mesh_, 
                species_, 
                i, 
                E_
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
    )
{
    constructModels();

    particleFlux_.setSize(species.nSpecies());

    for (label i = 0; i < species.nSpecies(); ++i)
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

// void plasmaTransport::correct()
// {
//     const label nSpecies = species_.nSpecies();

//     for (label i = 0; i < nSpecies; ++i)
//     {
//         surfaceScalarField phiE = fvc::flux(E_);
//         transportModels_[i].correct(phiE);

//         // Build the matrix
//         tmp<fvScalarMatrix> tEqn = transportModels_[i].nEqn();
//         fvScalarMatrix& nEqn = tEqn.ref();

//         // --- Lambda for the Debug Print ---
//         auto debugBreakdown = [&](const string& labelTag)
//         {
//             label eIdx = species_.speciesID("e");
//             volScalarField& n = species_.numberDensity(eIdx);
//             const scalarField& nCurr = n.internalField();
//             const scalarField& diag = nEqn.diag();
//             const scalarField& upper = nEqn.upper();
//             const scalarField& lower = nEqn.lower();
//             const scalarField& source = nEqn.source();
//             const labelList& owner = n.mesh().owner();
//             const labelList& neighbour = n.mesh().neighbour();

//             Pout << "--- [" << labelTag << "] Breakdown Species [" << i << "] ---" << endl;

//             for (label celli = 0; celli <= 5; ++celli) // Small range for readability
//             {
//                 scalar leftContrib = 0;
//                 scalar rightContrib = 0;

//                 for (label facei : n.mesh().cells()[celli])
//                 {
//                     if (n.mesh().isInternalFace(facei))
//                     {
//                         if (owner[facei] == celli) 
//                             rightContrib += upper[facei] * nCurr[neighbour[facei]];
//                         else 
//                             leftContrib += lower[facei] * nCurr[owner[facei]];
//                     }
//                 }

//                 scalar diagTerm = diag[celli] * nCurr[celli];
//                 scalar residual = (diagTerm + leftContrib + rightContrib) - source[celli];

//                 Pout << "Cell " << celli << " | Diag: " << diagTerm 
//                      << " | L: " << leftContrib << " | R: " << rightContrib 
//                      << " | Src: " << source[celli] 
//                      << " | Res: " << residual << endl;
//             }
//         };

//         // 1. Check BEFORE solve (shows the current imbalance/error)
//         debugBreakdown("BEFORE SOLVE");

//         // 2. Solve
//         nEqn.solve();

//         // 3. Check AFTER solve (residual should be ~0)
//         debugBreakdown("AFTER SOLVE");
//     }
// }

void plasmaTransport::correct(const bool firstIter, const bool finalIter)
{
    const label nSpecies = species_.nSpecies();

    label eIdx = species_.speciesID("e");
    label iIdx = species_.speciesID("pIon");
    
    // 1. Electric field flux calculation
    surfaceScalarField phiE = -fvc::snGrad(ePotential_) * mesh_.magSf();

    // 2. Update transport properties (Mobility/Diffusivity)
    for (label i = 0; i < nSpecies; ++i)
    {
        transportModels_[i].correct(phiE);
    }

    // 3. Calculate Electric field magnitude
    volScalarField Emag(mag(E_));

    // --- SAFETY BLOCK ---
    // Use a field-aware max to protect internal cells AND processor patches.
    // This prevents SIGFPE division by zero after mesh redistribution.
    volScalarField safeEmag = max(Emag, dimensionedScalar("minE", Emag.dimensions(), 1.0));

    // 4. Ionization coefficient (alpha) calculation
    dimensionedScalar E_const("E_const", Emag.dimensions(), 2.73e7);
    dimensionedScalar E_pow_const("E_pow_const", pow(Emag.dimensions(), 3), 4.3666e26);
    dimensionedScalar alpha_scale("alpha_scale", dimensionSet(0, -1, 0, 0, 0, 0, 0), 1.0);

    volScalarField alpha = 
        (1.1944e6 + E_pow_const / pow(safeEmag, 3)) 
        * exp(-E_const / safeEmag)
        * alpha_scale;
    
    alpha.correctBoundaryConditions();

    // 5. Source term assembly
    dimensionedScalar eta("eta", alpha.dimensions(), 340.75);
    volScalarField alphaEff = alpha - eta;
    volScalarField veMag = mag(transportModels_[eIdx].driftVelocity());

    k_eff_ = alphaEff * veMag;
    k_eff_.correctBoundaryConditions();

    volScalarField explicitSource = k_eff_ * species_.numberDensity(eIdx);

    // 6. Solve Continuity Equations
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

    for (label i = 0; i < nSpecies; ++i)
    {
        transportModels_[i].updateParticleFlux(particleFlux_[i]);
    }
}

void plasmaTransport::correctModels()
{
    const label nSpecies = species_.nSpecies();


    surfaceScalarField phiE = -fvc::snGrad(ePotential_) * mesh_.magSf();

    // Update transport properties
    for (label i = 0; i < nSpecies; ++i)
    {
        transportModels_[i].correct(phiE);
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












