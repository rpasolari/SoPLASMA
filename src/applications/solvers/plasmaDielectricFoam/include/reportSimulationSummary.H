/*---------------------------------------------------------------------------*\
  File: reportSimulationSummary.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Summarizes the complete simulation configuration to the console
      before the time loop begins.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

Info<<
"=============================================================================="
<< nl << "                      PLASMA SIMULATION CONFIGURATION" << nl <<
"==============================================================================" 
<< nl << endl;

// Time Stepping Controls
{
    Info<< "Time Stepping Control:" << nl;
    
    if (timeControl.adjustTimeStep())
    {
        Info<< "    Strategy:       ADJUSTABLE (Computed based on constraints)"
        << nl;
        
        if (timeControl.limitDielectricRelaxationRatio())
        {
            Info<< "    Constraint:     Dielectric Relaxation Ratio < " 
                << timeControl.maxDielectricRelaxationRatio() << nl;
        }

        if (timeControl.limitSpeciesCo())
        {
            Info<< "    Constraint:     Species Courant Number      < " 
                << timeControl.maxSpeciesCo() 
                << " (Species: " << timeControl.speciesName() << ")" << nl;
        }
    }
    else
    {
        Info<< "    Strategy:       CONSTANT" << nl;
        Info<< "    Delta T:        " << runTime.deltaT().value() << " s" << nl;
    }
    Info<< endl;
}

// // Adaptive Mesh Refinement (AMR) Configuration
// {
//     Info<< "Adaptive Mesh Refinement (AMR):" << nl;

//     if (enableAMR)
//     {
//         Info<< "    Status:         ENABLED" << nl;
//         Info<< "    Refine Interval: " << refineInterval << " steps" << nl;
//         Info<< "    Unrefine Interval: " << unrefineInterval << " steps" << nl;
//         Info<< "    Balance Interval: " << balanceInterval << " steps" << nl;

//         if (restrictAMRToPlasma)
//         {
//             Info<< "    Mode:           RESTRICTED (Plasma region only)" << nl;
//         }
//         else
//         {
//             Info<< "    Mode:           MULTI-REGION" << nl;
            
//             word syncStr = conformalAMRSync ? "ENABLED (Conformal)" : "DISABLED";
//             Info<< "    Sync Strategy:  " << syncStr << nl;

//             // Report which regions are participating in the sync
//             if (conformalAMRSync)
//             {
//                 string syncRegions = gasMesh().name();
//                 forAll(dielectricRegions, i)
//                 {
//                     if (touchesGas[i])
//                     {
//                         syncRegions += ", " + dielectricRegions[i].name();
//                     }
//                 }
//                 Info<< "    Sync Regions:   [" << syncRegions << "]" << nl;
//             }
//         }
//     }
//     else
//     {
//         Info<< "    Status:         DISABLED" << nl;
//     }
//     Info<< endl;
// }

// Electrostatics (Poisson & Region Coupling)
{
    Info<< "Electrostatics (Poisson & Region Coupling):" << nl;

    if (coupled)
    {
        Info<< "    Mode:           MONOLITHIC (Implicit coupling)" << nl;
        Info<< "    Non-Ortho Corr: " << nonOrthoCorrCoupled 
            << " (Global Max)" << nl;
    }
    else if (nTotalRegions > 1)
    {
        Info<< "    Mode:           SEGREGATED (Iterative multi-region)" << nl;
        Info<< "    Max Outer Iter: " << ePotentialNonCoupledLoopMaxIter << nl;
        Info<< "    Tolerance:      " << ePotentialNonCoupledTolerance 
            << " (Gas Region)" << nl;
    }
    else
    {
        // Only 1 region exists
        Info<< "    Mode:           SINGLE-REGION (Stand-alone)" << nl;
        Info<< "    Non-Ortho Corr: " << nNonOrthoCorr << nl;
    }

    Info<< "    Poisson Solver: " << poissonSolver << nl << endl;
}

// Species Transport Physics and Numerics
{
    Info<< "Species Transport Models:" << nl << nl;
    
    auto pad = [](const string& str, const label width) -> string
    {
        string s = str;
        if (s.size() >= width) 
        {
            return s.substr(0, width-1) + " "; 
        }
        return s + string(width - s.size(), ' ');
    };

    Info<< pad("Species", 15).c_str() 
        << pad("Charge", 10).c_str() 
        << pad("Transport Model", 20).c_str() 
        << pad("Scheme", 25).c_str() 
        << pad("Advection", 18).c_str()
        << pad("Mobility", 15).c_str()
        << pad("Diffusivity", 15).c_str()
        << nl;

    Info<< pad("-------", 15).c_str() 
        << pad("------", 10).c_str() 
        << pad("---------------", 20).c_str() 
        << pad("------", 25).c_str() 
        << pad("---------", 18).c_str()
        << pad("--------", 15).c_str()
        << pad("-----------", 15).c_str()
        << nl;

    for (label i = 0; i < species.nSpecies(); ++i)
    {
        const word& sName = species.speciesNames()[i];
        scalar q = species.speciesChargeNumber(i);
        const dictionary& sDict = species.speciesDict(i);

        word modelType =
            sDict.lookupOrDefault<word>("transportModel", "unknown");
        word mobModel = sDict.lookupOrDefault<word>("mobilityModel", "none");
        word diffModel =
            sDict.lookupOrDefault<word>("diffusivityModel", "none");
        
        word schemeStr = "n/a";
        word advStr = "implicit"; 

        if (modelType == "driftDiffusion")
        {
            schemeStr =
            sDict.lookupOrDefault<word>("driftDiffusionFluxScheme", "standard");
            
            if (schemeStr == "ScharfetterGummel")
            {
                advStr = "implicit (SG)";
            }
            else
            {
                if (sDict.found("advectionMode"))
                {
                    word mode;
                    sDict.lookup("advectionMode") >> mode;
                    
                    if (mode == "explicit") 
                        advStr = "explicit";
                }
            }
        }
        else
        {
            schemeStr = "standard"; 
            if (sDict.found("advectionMode"))
            {
                 word mode;
                 sDict.lookup("advectionMode") >> mode;
                 
                 if (mode == "explicit") 
                     advStr = "explicit";
            }
        }

        Info<< pad(sName, 15).c_str() 
            << pad(Foam::name(q), 10).c_str() 
            << pad(modelType, 20).c_str() 
            << pad(schemeStr, 25).c_str() 
            << pad(advStr, 18).c_str()
            << pad(mobModel, 15).c_str()
            << pad(diffModel, 15).c_str()
            << nl;
    }
    Info<< endl;
}

Info<<
"=============================================================================="
<< nl << "                      END OF CONFIGURATION REPORT" << nl <<
"==============================================================================" 
<< endl;
