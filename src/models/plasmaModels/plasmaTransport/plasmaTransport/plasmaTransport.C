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
#include "plasmaProfiler.H"

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

void plasmaTransport::correct()
{
    // Update transport coefficients for all models
    for (label i = 0; i < species_.nSpecies(); ++i)
    {
        transportModels_[i].correct();
    }
}

// void plasmaTransport::correct(const bool firstIter, const bool finalIter)
// {
//     const label eIdx = species_.speciesID("e");
//     const label iIdx = species_.speciesID("pIon");
    
//     const volScalarField& Emag = species_.Emag();

//     volScalarField safeEmag = max(Emag, dimensionedScalar("minE", Emag.dimensions(), 1e-6));
    
//     // 3. Ionization coefficient (alpha) calculation for 0.1 torr
//     // A*p = 0.1 * 34 = 3.4 cm^-1 = 340 m^-1 [cite: 242, 243]
//     dimensionedScalar Ap("Ap", dimensionSet(0, -1, 0, 0, 0, 0, 0), 340.0);
    
//     // Stripping dimensions for math functions
//     dimensionedScalar E_ref("E_ref", safeEmag.dimensions(), 1.0);
//     volScalarField E_num = safeEmag / E_ref;

//     // Calculation: Ap is m^-1, exp part is dimensionless. Result is m^-1.
//     // The factor 0.1 converts E_num (V/m) to the paper's units (V/cm-torr)
//     volScalarField alpha = Ap * exp(-16.0 / pow(0.1 * E_num, 0.4));
//     alpha.correctBoundaryConditions();

//     // 4. Electron drift velocity calculation
//     const driftDiffusion& eModel = refCast<const driftDiffusion>(transportModels_[eIdx]);
//     const auto& mobE = eModel.mobility();

//     volScalarField veMag
//     (
//         IOobject("veMag", mesh_.time().timeName(), mesh_),
//         mesh_,
//         dimensionedScalar("zeroVe", dimensionSet(0, 1, -1, 0, 0, 0, 0), 0.0)
//     );

//     if (mobE.isUniform())
//     {
//         veMag = mobE.muValue() * Emag;
//     }
//     else
//     {
//         veMag = mobE.mu() * Emag;
//     }

//     // 5. Source term assembly
//     k_eff_ = alpha * veMag;
//     k_eff_.correctBoundaryConditions();
//     volScalarField explicitSource = k_eff_ * species_.numberDensity(eIdx);

//     // 6. Solve Continuity Equations
//     tmp<fvScalarMatrix> tEqne = transportModels_[eIdx].nEqn();
//     tEqne.ref() -= explicitSource; 
//     tEqne.ref().solve();
    
//     tmp<fvScalarMatrix> tEqni = transportModels_[iIdx].nEqn();
//     tEqni.ref() -= explicitSource; 
//     tEqni.ref().solve();

//     // 7. Density Limiting
//     for (const label id : {eIdx, iIdx})
//     {
//         volScalarField& n = species_.numberDensity(id);
//         n = max(n, dimensionedScalar("nMin", n.dimensions(), 1e10));
//         n.correctBoundaryConditions();
//     }

//     for (const label i : species_.mobileSpeciesIDs())
//     {
//         transportModels_[i].updateParticleFlux(particleFlux_[i]);
//     }
// }
// void plasmaTransport::solve()
// {
//     const label nSpecies = species_.nSpecies();
//     const label eIdx = species_.speciesID("e");
//     const label iIdx = species_.speciesID("pIon");
    
//     // Electric field flux calculation
//     phiE_ = -fvc::snGrad(ePotential_) * mesh_.magSf();

//     // Calculate Electric field magnitude
//     volScalarField Emag(mag(E_));

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

//     volScalarField explicitSource = k_eff_ * species_.numberDensity(eIdx);

//     // Solve Continuity Equations
//     tmp<fvScalarMatrix> tEqne = transportModels_[eIdx].nEqn();
//     fvScalarMatrix& nEqne = tEqne.ref();
//     nEqne -= explicitSource; 
    
//     tmp<fvScalarMatrix> tEqni = transportModels_[iIdx].nEqn();
//     fvScalarMatrix& nEqni = tEqni.ref();
//     nEqni -= explicitSource; 

//     nEqni.solve();
//     nEqne.solve();

//     volScalarField& ne = species_.numberDensity(eIdx);
//     dimensionedScalar neMin("neMin", ne.dimensions(), 1e10);
//     ne = max(ne, neMin);
//     ne.correctBoundaryConditions();

//     volScalarField& ni = species_.numberDensity(iIdx);
//     dimensionedScalar niMin("niMin", ni.dimensions(), 1e10);
//     ni = max(ni, niMin);
//     ni.correctBoundaryConditions();

//     for (const label i : species_.mobileSpeciesIDs())
//     {
//         transportModels_[i].updateParticleFlux(particleFlux_[i]);
//     }
// }

void plasmaTransport::solve()
{
    //
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












