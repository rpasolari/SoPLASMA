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
                E_,
                ePotential_
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
    particleFlux_()
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
//         // Update coefficients (mobility, diffusivity, driftVelocity etc.)
//         transportModels_[i].correct();

//         // Build the matrix
//         tmp<fvScalarMatrix> tEqn = transportModels_[i].nEqn();
//         fvScalarMatrix& nEqn = tEqn.ref();

//         // Solve the transport equation
//         nEqn.solve();

//         // Update fluxes
//         particleFlux_[i] = transportModels_[i].particleFlux();
//         Gamma_[i] = fvc::reconstruct(particleFlux_[i]);
//     }
// }

void plasmaTransport::correct()
{
    const label nSpecies = species_.nSpecies();
    label eIdx = species_.speciesID("e");
    label iIdx = species_.speciesID("pIon");

    // Update transport properties
    for (label i = 0; i < nSpecies; ++i)
    {
        transportModels_[i].correct();
    }

    // 2. Extract Electric Field magnitude and normalize for empirical fits
    const volScalarField& ne = species_.numberDensity(eIdx);
    volScalarField Emag(mag(E_));
    
    dimensionedScalar E_unit("E_unit", Emag.dimensions(), 1.0);
    volScalarField E_num = max(Emag, E_unit) / E_unit;

    // 3. Calculate Ionization (alpha) and Attachment (eta) coefficients
    volScalarField alpha = (1.1944e6 + 4.3666e26/pow(E_num, 3)) * exp(-2.73e7/E_num) 
                         * dimensionedScalar("invM", dimensionSet(0,-1,0,0,0,0,0), 1.0);
    
    dimensionedScalar eta("eta", dimensionSet(0,-1,0,0,0,0,0), 340.75);

    // 4. Calculate Drift Velocity magnitude and effective rate
    volScalarField veMag = mag(transportModels_[eIdx].driftVelocity());
    volScalarField k_eff = (alpha - eta) * veMag;

    // 5. Solve Electron Continuity Equation
    {
        tmp<fvScalarMatrix> tEqne = transportModels_[eIdx].nEqn();
        fvScalarMatrix& nEqne = tEqne.ref();

        // Use SuSp for the source term S = k_eff * ne
        nEqne -= fvm::SuSp(k_eff, ne); 

        nEqne.solve();
    }

    // 6. Solve Positive Ion Equation (Immobile ions)
    {
        const volScalarField& neNew = species_.numberDensity(eIdx);
        
        tmp<fvScalarMatrix> tEqni = transportModels_[iIdx].nEqn();
        fvScalarMatrix& nEqni = tEqni.ref();

        // Ion production rate 
        nEqni -= k_eff * neNew;

        nEqni.solve();
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












