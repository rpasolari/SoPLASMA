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

// Separators and helpers

const std::string HL(72, '*');
const std::string sl(72, '-');

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

// Header

Info<< nl
    << HL.c_str() << nl
    << "##" << nl
    << "##    SoPLASMA  —  simulation configuration" << nl
    << "##" << nl
    << HL.c_str() << nl << endl;

// ---- Electromagnetics ----

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

// Electromagnetics: Dielectric region

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

// ---- Species & Transport ----

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
    + col("mu / D", 30);

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
        + col(muD,                     30);

    Info<< row.c_str() << nl;
}

// Footer

Info<< nl
    << HL.c_str() << nl
    << "##" << nl
    << "##    end of configuration" << nl
    << "##" << nl
    << HL.c_str() << nl << endl;

