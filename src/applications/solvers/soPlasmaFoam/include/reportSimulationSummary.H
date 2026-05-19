/*---------------------------------------------------------------------------*\
  File: reportSimulationSummary.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Summarizes the complete simulation configuration to the console
      before the time loop begins.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

// ----- SEPARATORS AND HELPERS -----

const std::string HL(72, '*');
const std::string sl(72, '-');

auto field = [](const std::string& lbl, const std::string& val, int w=28)
{
    std::string p = lbl + ":";
    while (static_cast<int>(p.size()) < w) p += ' ';
    return p + val;
};

auto fieldPair = [](
    const std::string& l1, const std::string& v1,
    const std::string& l2, const std::string& v2,
    int w1=18, int w2=22)
{
    std::string a = l1 + ":";
    while (static_cast<int>(a.size()) < w1) a += ' ';
    a += v1;
    while (static_cast<int>(a.size()) < w1 + w2) a += ' ';
    std::string b = l2 + ":";
    while (static_cast<int>(b.size()) < w1) b += ' ';
    return a + b + v2;
};

auto labelListStr = [](const labelList& ids) -> std::string
{
    if (ids.empty()) return "(none)";
    std::string s = "(";
    forAll(ids, i) { if (i) s += " "; s += std::to_string(ids[i]); }
    return s + ")";
};

// ----- HEADER -----

Info<< nl
    << HL.c_str() << nl
    << "##" << nl
    << "##    SoPLASMA  —  simulation configuration" << nl
    << "##" << nl
    << HL.c_str() << nl << endl;

// ----- ELECTROMAGNETICS : GAS REGION -----

Info<< "  ELECTROMAGNETICS — gas region" << nl
    << "  " << sl.c_str() << nl
    << "  " << fieldPair("Model",         em->type(),
                         "Poisson scheme", em->PoissonScheme())              << nl
    << "  " << fieldPair("E-field scheme", em->EScheme(),
                         "Non-ortho corr.", Foam::name(em->nNonOrthCorr())) << nl
    << "  " << fieldPair("epsilonR",       Foam::name(em->epsilonR()),
                         "epsilon [F/m]",  Foam::name(em->epsilon().value())) << nl;

const multiRegionPoissonModel* mrp = isA<multiRegionPoissonModel>(em());
if (mrp)
{
    const word coupling = mrp->coupled() ? "monolithic" : "segregated";
    Info<< "  " << fieldPair("Coupling",     coupling,
                             "Dielectrics",  Foam::name(mrp->nDielectrics())) << nl;
}
Info<< nl;

// ----- ELECTROMAGNETICS : DIELECTRIC REGIONS -----

if (mrp && mrp->nDielectrics() > 0)
{
    Info<< "  ELECTROMAGNETICS — dielectric regions" << nl
        << "  " << sl.c_str() << nl;

    forAll(mrp->dielectrics(), i)
    {
        const dielectricRegion& reg = mrp->dielectric(i);
        Info<< "  " << fieldPair(
                "[" + std::to_string(i) + "] " + reg.mesh().name(), "",
                "epsilonR", Foam::name(reg.epsilonR())) << nl;
    }
    Info<< nl;
}

// ----- SPECIES & TRANSPORT -----

auto col = [](const std::string& s, int w) -> std::string
{
    if (static_cast<int>(s.size()) >= w) return s.substr(0, w-1) + " ";
    return s + std::string(w - s.size(), ' ');
};

Info<< "  SPECIES & TRANSPORT" << nl
    << "  " << sl.c_str() << nl
    << "  " << fieldPair("Background gas",  species.backgroundName(),
                         "N [1/m3]",
                         Foam::name(species.backgroundDensity().value()))  << nl
    << "  " << fieldPair("Active species",
                         Foam::name(species.nSpecies()),
                         "Electron ID",
                         Foam::name(species.electronSpeciesID()))          << nl
    << "  " << fieldPair("Mobile IDs",
                         labelListStr(species.mobileSpeciesIDs()),
                         "Immobile IDs",
                         labelListStr(species.immobileSpeciesIDs()))       << nl
    << nl;

// ----- species table -----

const std::string hdr =
    "  "
    + col("ID",        4)
    + col("Name",      8)
    + col("Role",      13)
    + col("Z",         5)
    + col("Mass [kg]", 14)
    + col("minN [1/m3]", 13)
    + col("Transport", 17)
    + col("Flux scheme", 20)
    + col("mu / D", 20);

Info<< hdr.c_str() << nl
    << "  " << sl.c_str() << nl;

forAll(species.speciesNames(), i)
{
    const scalar Z = species.speciesChargeNumber(i);
    const scalar m = species.speciesMass(i).value();
    const scalar minN = species.speciesMinNumberDensity(i);


    std::string role;
    if      (i == species.electronSpeciesID()) role = "electron";
    else if (Z > 0)                            role = "pos. ion";
    else if (Z < 0)                            role = "neg. ion";
    else                                       role = "neutral";

    const plasmaTransportModel& tm = transport.model(i);

    std::string scheme  = "-";
    std::string muD     = "-";

    const driftDiffusion* dd = isA<driftDiffusion>(tm);
    if (dd)
    {
        scheme = dd->fluxScheme();
        muD    = dd->mobility().modelType()
               + " / "
               + dd->diffusivity().modelType();
    }

    const std::string row =
        "  "
        + col(std::to_string(i),       4)
        + col(species.speciesName(i),  8)
        + col(role,                    13)
        + col(Foam::name(Z),           5)
        + col(Foam::name(m),           14)
        + col(Foam::name(minN),        13)
        + col(tm.type(),               17)
        + col(scheme,                  20)
        + col(muD,                     20);

    Info<< row.c_str() << nl;
}

// ----- FOOTER -----

Info<< nl
    << HL.c_str() << nl
    << "##" << nl
    << "##    end of configuration" << nl
    << "##" << nl
    << HL.c_str() << nl << endl;

// Info<<
// "=============================================================================="
// << nl << "                      PLASMA SIMULATION CONFIGURATION" << nl <<
// "==============================================================================" 
// << nl << endl;

// // Time Stepping Controls
// {
//     Info<< "Time Stepping Control:" << nl;
    
//     if (timeControl.adjustTimeStep())
//     {
//         Info<< "    Strategy:       ADJUSTABLE (Computed based on constraints)"
//         << nl;
        
//         if (timeControl.limitDielectricRelaxationRatio())
//         {
//             Info<< "    Constraint:     Dielectric Relaxation Ratio < " 
//                 << timeControl.maxDielectricRelaxationRatio() << nl;
//         }

//         if (timeControl.limitSpeciesCo())
//         {
//             Info<< "    Constraint:     Species Courant Number      < " 
//                 << timeControl.maxSpeciesCo() 
//                 << " (Species: " << timeControl.speciesName() << ")" << nl;
//         }
//     }
//     else
//     {
//         Info<< "    Strategy:       CONSTANT" << nl;
//         Info<< "    Delta T:        " << runTime.deltaT().value() << " s" << nl;
//     }
//     Info<< endl;
// }

// // // Adaptive Mesh Refinement (AMR) Configuration
// // {
// //     Info<< "Adaptive Mesh Refinement (AMR):" << nl;

// //     if (enableAMR)
// //     {
// //         Info<< "    Status:         ENABLED" << nl;
// //         Info<< "    Refine Interval: " << refineInterval << " steps" << nl;
// //         Info<< "    Unrefine Interval: " << unrefineInterval << " steps" << nl;
// //         Info<< "    Balance Interval: " << balanceInterval << " steps" << nl;

// //         if (restrictAMRToPlasma)
// //         {
// //             Info<< "    Mode:           RESTRICTED (Plasma region only)" << nl;
// //         }
// //         else
// //         {
// //             Info<< "    Mode:           MULTI-REGION" << nl;
            
// //             word syncStr = conformalAMRSync ? "ENABLED (Conformal)" : "DISABLED";
// //             Info<< "    Sync Strategy:  " << syncStr << nl;

// //             // Report which regions are participating in the sync
// //             if (conformalAMRSync)
// //             {
// //                 string syncRegions = gasMesh().name();
// //                 forAll(dielectricRegions, i)
// //                 {
// //                     if (touchesGas[i])
// //                     {
// //                         syncRegions += ", " + dielectricRegions[i].name();
// //                     }
// //                 }
// //                 Info<< "    Sync Regions:   [" << syncRegions << "]" << nl;
// //             }
// //         }
// //     }
// //     else
// //     {
// //         Info<< "    Status:         DISABLED" << nl;
// //     }
// //     Info<< endl;
// // }

// // Electrostatics (Poisson & Region Coupling)
// {
//     Info<< "Electrostatics (Poisson & Region Coupling):" << nl;

//     if (coupled)
//     {
//         Info<< "    Mode:           MONOLITHIC (Implicit coupling)" << nl;
//         Info<< "    Non-Ortho Corr: " << nonOrthoCorrCoupled 
//             << " (Global Max)" << nl;
//     }
//     else if (nTotalRegions > 1)
//     {
//         Info<< "    Mode:           SEGREGATED (Iterative multi-region)" << nl;
//         Info<< "    Max Outer Iter: " << ePotentialNonCoupledLoopMaxIter << nl;
//         Info<< "    Tolerance:      " << ePotentialNonCoupledTolerance 
//             << " (Gas Region)" << nl;
//     }
//     else
//     {
//         // Only 1 region exists
//         Info<< "    Mode:           SINGLE-REGION (Stand-alone)" << nl;
//         Info<< "    Non-Ortho Corr: " << nNonOrthoCorr << nl;
//     }

//     Info<< "    Poisson Solver: " << poissonSolver << nl << endl;
// }

// // Species Transport Physics and Numerics
// {
//     Info<< "Species Transport Models:" << nl << nl;
    
//     auto pad = [](const string& str, const label width) -> string
//     {
//         string s = str;
//         if (s.size() >= width) 
//         {
//             return s.substr(0, width-1) + " "; 
//         }
//         return s + string(width - s.size(), ' ');
//     };

//     Info<< pad("Species", 15).c_str() 
//         << pad("Charge", 10).c_str() 
//         << pad("Transport Model", 20).c_str() 
//         << pad("Scheme", 25).c_str() 
//         << pad("Advection", 18).c_str()
//         << pad("Mobility", 15).c_str()
//         << pad("Diffusivity", 15).c_str()
//         << nl;

//     Info<< pad("-------", 15).c_str() 
//         << pad("------", 10).c_str() 
//         << pad("---------------", 20).c_str() 
//         << pad("------", 25).c_str() 
//         << pad("---------", 18).c_str()
//         << pad("--------", 15).c_str()
//         << pad("-----------", 15).c_str()
//         << nl;

//     for (label i = 0; i < species.nSpecies(); ++i)
//     {
//         const word& sName = species.speciesNames()[i];
//         scalar q = species.speciesChargeNumber(i);
//         const dictionary& sDict = species.speciesDict(i);

//         word modelType =
//             sDict.lookupOrDefault<word>("transportModel", "unknown");
//         word mobModel = sDict.lookupOrDefault<word>("mobilityModel", "none");
//         word diffModel =
//             sDict.lookupOrDefault<word>("diffusivityModel", "none");
        
//         word schemeStr = "n/a";
//         word advStr = "implicit"; 

//         if (modelType == "driftDiffusion")
//         {
//             schemeStr =
//             sDict.lookupOrDefault<word>("driftDiffusionFluxScheme", "standard");
            
//             if (schemeStr == "ScharfetterGummel")
//             {
//                 advStr = "implicit (SG)";
//             }
//             else
//             {
//                 if (sDict.found("advectionMode"))
//                 {
//                     word mode;
//                     sDict.lookup("advectionMode") >> mode;
                    
//                     if (mode == "explicit") 
//                         advStr = "explicit";
//                 }
//             }
//         }
//         else
//         {
//             schemeStr = "standard"; 
//             if (sDict.found("advectionMode"))
//             {
//                  word mode;
//                  sDict.lookup("advectionMode") >> mode;
                 
//                  if (mode == "explicit") 
//                      advStr = "explicit";
//             }
//         }

//         Info<< pad(sName, 15).c_str() 
//             << pad(Foam::name(q), 10).c_str() 
//             << pad(modelType, 20).c_str() 
//             << pad(schemeStr, 25).c_str() 
//             << pad(advStr, 18).c_str()
//             << pad(mobModel, 15).c_str()
//             << pad(diffModel, 15).c_str()
//             << nl;
//     }
//     Info<< endl;
// }

// Info<<
// "=============================================================================="
// << nl << "                      END OF CONFIGURATION REPORT" << nl <<
// "==============================================================================" 
// << endl;
