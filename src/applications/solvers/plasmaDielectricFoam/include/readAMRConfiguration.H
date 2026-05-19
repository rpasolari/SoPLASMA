
/*---------------------------------------------------------------------------*\
  File: readAMRConfiguration.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

Description:
    Performs configuration checks and initializations for multi-region Adaptive 
    Mesh Refinement (AMR). This file validates solver controls and identifies 
    coupling interfaces to ensure safety during implicit coupling.

    The multi-region refinement logic implemented in this toolkit is based on 
    the methodology developed in the blastAMR project:
    https://github.com/STFS-TUDa/blastAMR

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

const dictionary& AMRDict = plasmaSimulationControlsDict.subDict("AMR");

Switch enableAMR = AMRDict.lookupOrDefault<bool>("enable", false);
Switch restrictAMRToPlasma = false;
Switch conformalAMRSync = false;

adaptiveFvMesh& gasAMRMesh = refCast<adaptiveFvMesh>(gasMesh());

label totalRegions = dielectricRegions.size() + 1;
autoPtr<volScalarField> refineFieldGasPtr;
PtrList<volScalarField> refineFieldDielectricList(dielectricRegions.size());
boolList touchesGas(dielectricRegions.size(), false);

label refineInterval = 1;
label unrefineInterval = 1;
label balanceInterval = 1;
if (enableAMR)
{
    IOdictionary dynamicMeshDict
    (
        IOobject
        (
            "dynamicMeshDict",
            gasMesh().time().constant(),
            gasMesh(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    refineInterval =
        dynamicMeshDict.lookupOrDefault<label>("refineInterval", 1);
    unrefineInterval = 
        dynamicMeshDict.lookupOrDefault<label>("unrefineInterval", 1);
    balanceInterval = 
       dynamicMeshDict.lookupOrDefault<label>("balanceInterval", refineInterval);
}

// AMR Controls Reading
if (enableAMR)
{
    if (AMRDict.found("AMRControls"))
    {
        const dictionary& amrControls = AMRDict.subDict("AMRControls");

        restrictAMRToPlasma = amrControls.lookupOrDefault<bool>
        (
            "restrictAMRToPlasma",
            false
        );

        conformalAMRSync = amrControls.lookupOrDefault<bool>
        (
            "conformalAMRSync",
            false
        );
    }

    Info<< "\nAdaptive Mesh Refinement (AMR) is ENABLED." << endl;

    if (restrictAMRToPlasma)
    {
        IOWarningInFunction(plasmaSimulationControlsDict)
            << "AMR is restricted to the MASTER (GAS) region only." << endl
            << "1. Ensure 'dynamicMeshDict' is only active in the gas region."
            << endl
            << "2. PROTECTION: Use 'protectedPatches' in dynamicMeshDict." 
            << endl
            << "3. RISK: Non-conformal interfaces crash implicit coupling."
            << endl
            << "See: SoPLASMA/docs/simulationManual/AMR.md" 
            << endl;
    }
    else
    {        
        if (conformalAMRSync)
        {
            Info<< "1. 'conformalAMRSync' is ENABLED." << endl;
            Info<< "2. CRITICAL: 'dynamicMeshDict' settings (maxRefinement, "
                << "nBufferLayers) MUST be identical across regions." << endl;
        }
        else
        {
            IOWarningInFunction(plasmaSimulationControlsDict)
                << "'conformalAMRSync' is DISABLED." << endl
                << "1. Implicitly coupled regions (mappedWall) at HIGH RISK." 
                << endl
                << "2. Enable 'conformalAMRSync' or protect interface patches."
                << endl
                << "See: SoPLASMA/docs/simulationManual/AMR.md" 
                << endl;
        }
    }

    if (Pstream::parRun())
    {
        Info<< "\nPARALLEL RUN DETECTED:" << endl
            << "-> Use 'TAKA' in decomposeParDict for interfaces."
            << endl
            << "See: SoPLASMA/docs/simulationManual/AMR.md" << endl;
    }
}

if (enableAMR && !restrictAMRToPlasma && totalRegions > 1 && conformalAMRSync)
{
    // Initialize Gas refineField that will be passed in neighbours
    refineFieldGasPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "refineField",
                runTime.timeName(),
                gasMesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            gasMesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
    
    forAll(dielectricRegions, i)
    {
        const fvBoundaryMesh& boundary = dielectricRegions[i].boundary();
        touchesGas[i] = false; 

        forAll(boundary, patchi)
        {
            const polyPatch& pp = boundary[patchi].patch();

            if (isA<mappedPatchBase>(pp))
            {
                const mappedPatchBase& mpp = refCast<const mappedPatchBase>(pp);

                if (mpp.sampleRegion() == gasMesh().name())
                {
                    touchesGas[i] = true;
                    break; 
                }
            }
        }

        if (touchesGas[i])
        {
            refineFieldDielectricList.set
            (
                i,
                new volScalarField
                (
                    IOobject
                    (
                        "refineField",
                        runTime.timeName(),
                        dielectricRegions[i],
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    dielectricRegions[i],
                    dimensionedScalar("zero", dimless, 0.0)
                )
            );
        }
    }
}
