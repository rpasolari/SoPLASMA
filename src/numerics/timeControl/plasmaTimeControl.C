/*---------------------------------------------------------------------------*\
  File: plasmaTimeControl.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaTimeControl.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "plasmaTimeControl.H"
#include "plasmaTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

plasmaTimeControl::plasmaTimeControl(Time& runTime, const fvMesh& mesh)
:
    runTime_(runTime),
    mesh_(mesh),
    dict_(dictionary::null),
    adjustTimeStep_(false),
    maxDeltaT_(GREAT),
    limitDielectricRelaxationRatio_(false),
    printDielectricRelaxationRatio_(false),
    maxDielectricRelaxationRatio_(1.0),
    limitSpeciesCo_(false),
    printSpeciesCo_(false),
    maxSpeciesConvectiveCo_(1.0),
    maxSpeciesDiffusiveCo_(1.0),
    courantSpeciesName_("e"),
    limitChemistryCo_(false),
    printChemistryCo_(false),
    maxChemistryCo_(1.0),
    limitVoltageRiseRate_(false),
    printVoltageRiseRate_(false),
    maxVoltageRiseRate_(GREAT),
    voltagePatchName_(""),
    prevPatchVoltage_(0.0)
{
    read();

    Info << "Plasma time control succesfully constructed." << nl << endl;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void plasmaTimeControl::read()
{
    IOdictionary plasmaDict
        (
            IOobject
            (
                "plasmaSimulationControls",
                runTime_.system(),
                runTime_,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );


    if (plasmaDict.found("plasmaTimeControl"))
    {
        dict_ = plasmaDict.subDict("plasmaTimeControl");

        adjustTimeStep_ =
            dict_.lookupOrDefault<Switch>("adjustTimeStep", false);

        maxDeltaT_ =
            dict_.lookupOrDefault<scalar>("maxDeltaT", GREAT);

        // Dielectric relaxation
        limitDielectricRelaxationRatio_ =
         dict_.lookupOrDefault<Switch>("limitDielectricRelaxationRatio", false);

        printDielectricRelaxationRatio_ =
         dict_.lookupOrDefault<Switch>("printDielectricRelaxationRatio", false);

        maxDielectricRelaxationRatio_ =
            dict_.lookupOrDefault<scalar>("maxDielectricRelaxationRatio", 1.0);

        // Species Courant
        limitSpeciesCo_ = 
            dict_.lookupOrDefault<Switch>("limitSpeciesCo", false);

        printSpeciesCo_ = 
            dict_.lookupOrDefault<Switch>("printSpeciesCo", false);

        if (limitSpeciesCo_ || printSpeciesCo_)
        {
            maxSpeciesConvectiveCo_ =
                dict_.lookupOrDefault<scalar>("maxSpeciesConvectiveCo", 1.0);

            maxSpeciesDiffusiveCo_ =
                dict_.lookupOrDefault<scalar>("maxSpeciesDiffusiveCo", 1.0);

            courantSpeciesName_ =
                dict_.lookupOrDefault<word>("courantSpeciesName", "e");
        }

        // Chemistry Courant
        limitChemistryCo_ = 
            dict_.lookupOrDefault<Switch>("limitChemistryCo", false);

        printChemistryCo_ = 
            dict_.lookupOrDefault<Switch>("printChemistryCo", false);

        maxChemistryCo_ = 
            dict_.lookupOrDefault<scalar>("maxChemistryCo", 1.0);

        // Voltage rise rate
        limitVoltageRiseRate_ =
            dict_.lookupOrDefault<Switch>("limitVoltageRiseRate", false);

        printVoltageRiseRate_ =
            dict_.lookupOrDefault<Switch>("printVoltageRiseRate", false);

        if (limitVoltageRiseRate_ || printVoltageRiseRate_)
        {
            maxVoltageRiseRate_ =
                dict_.lookupOrDefault<scalar>("maxVoltageRiseRate", 1e5);

            voltagePatchName_ =
                dict_.lookupOrDefault<word>("voltagePatchName", "");

            if (voltagePatchName_.empty())
            {
                FatalIOErrorInFunction(dict_)
                    << "'voltagePatchName' must be specified when "
                    << "'limitVoltageRiseRate' or 'printVoltageRiseRate' "
                    << "is true."
                    << nl << exit(FatalIOError);
            }
        }
    }
}

scalar plasmaTimeControl::patchVoltageAvg
(
    const plasmaTransport& transport
) const
{
    const volScalarField& phi = transport.species().em().ePotential();
    const label patchi = mesh_.boundaryMesh().findPatchID(voltagePatchName_);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch '" << voltagePatchName_ << "' not found in mesh."
            << nl << exit(FatalError);
    }

    const scalarField& phiPatch = phi.boundaryField()[patchi];
    return phiPatch.empty() ? 0.0 : gAverage(phiPatch);
}

void plasmaTimeControl::adjustDeltaT(const plasmaTransport& transport)
{
    read();

    const scalar currentDeltaT = runTime_.deltaTValue();
    const scalar eps0 = constant::plasma::epsilon0.value();
    scalar newDeltaT = maxDeltaT_;

    scalar maxSigma = 0.0;
    scalar maxConvFluxRate  = 0.0;
    scalar maxDiffFluxRate  = 0.0;
    scalar meanConvFluxRate  = 0.0;
    scalar meanDiffFluxRate  = 0.0;
    scalar maxKeff = 0.0;
    scalar voltageRiseRate  = 0.0;

    //  Dielectric relaxation (tau = epsilon / sigma)
    if (limitDielectricRelaxationRatio_ || printDielectricRelaxationRatio_)
    {
        tmp<volScalarField> tSigma = transport.electricalConductivity();
        const volScalarField& sigma = tSigma();
        
        maxSigma = gMax(mag(sigma)().primitiveField());
        scalar dielectricLimit = 
            (maxDielectricRelaxationRatio_ * eps0) / (maxSigma + VSMALL);

        if (limitDielectricRelaxationRatio_)
        {
            newDeltaT = min(newDeltaT, dielectricLimit);
        }
        tSigma.clear();
    }
    
    //  Species Courant Limit (convective and diffusive)
    if (limitSpeciesCo_ || printSpeciesCo_)
    {
        label sIdx = transport.species().speciesID(courantSpeciesName_);
        const volScalarField& n = transport.species().numberDensity(sIdx);
        const scalarField& nField = n.primitiveField();
        const scalarField& V = mesh_.V().field();

        // Convective
        const surfaceScalarField& convFlux = transport.convectiveFlux(sIdx);
        const scalarField sumConv
        (
            fvc::surfaceSum(mag(convFlux))().primitiveField()
        );
        const scalarField convRate
        (
            0.5 * sumConv / (nField * V + VSMALL)
        );
        maxConvFluxRate = gMax(convRate);
        if (limitSpeciesCo_)
        {
            newDeltaT = min
            (
                newDeltaT,
                maxSpeciesConvectiveCo_ / (maxConvFluxRate + VSMALL)
            );
        }

        // Diffusive
        const surfaceScalarField& diffFlux = transport.diffusiveFlux(sIdx);
        const scalarField sumDiff
        (
            fvc::surfaceSum(mag(diffFlux))().primitiveField()
        );
        const scalarField diffRate
        (
            0.5 * sumDiff / (nField * V + VSMALL)
        );
        maxDiffFluxRate = gMax(diffRate);

        if (limitSpeciesCo_)
        {
            newDeltaT = min
            (
                newDeltaT,
                maxSpeciesDiffusiveCo_ / (maxDiffFluxRate + VSMALL)
            );
        }
    }

    // Chemistry Limit
    if (limitChemistryCo_ || printChemistryCo_)
    {
        // const volScalarField& keff = transport.k_eff();
        // maxKeff = gMax(mag(keff)().primitiveField());
        
        // scalar chemistryLimit = maxChemistryCo_ / (maxKeff + VSMALL);

        // if (limitChemistryCo_)
        // {
        //     newDeltaT = min
        //     (
        //         newDeltaT,
        //         maxChemistryCo_ / (maxKeff + VSMALL)
        //     );
        // }
    }

    // Voltage rise rate
    if (limitVoltageRiseRate_ || printVoltageRiseRate_)
    {
        const scalar currentVoltage = patchVoltageAvg(transport);
        const scalar dV = mag(currentVoltage - prevPatchVoltage_);
        voltageRiseRate = dV / (runTime_.deltaT0Value() + VSMALL);

        if (limitVoltageRiseRate_)
        {
            const scalar dtLimit =
                maxVoltageRiseRate_ / (voltageRiseRate + VSMALL);

            newDeltaT = min(newDeltaT, dtLimit);
        }

        prevPatchVoltage_ = currentVoltage;
    }

    // Apply and clamp
    if (adjustTimeStep_)
    {
        reduce(newDeltaT, minOp<scalar>());

        if (newDeltaT > currentDeltaT)
        {
            newDeltaT = min(newDeltaT, currentDeltaT * 1.2);
        }

        runTime_.setDeltaT(newDeltaT);
    }

    // Report
    const scalar actualDeltaT = runTime_.deltaTValue();

    word   bindingName  = "maxDeltaT";

    auto fmtLine = [](
        const std::string& label,
        scalar value,
        bool   limited,
        scalar limitVal,
        bool   isBinding,
        int    lw = 26) -> std::string
    {
        std::string line = "  " + label + ":";
        while (static_cast<int>(line.size()) < lw) line += ' ';
        line += Foam::name(value).c_str();

        if (limited)
        {
            while (static_cast<int>(line.size()) < lw + 14) line += ' ';
            line += "[ max: " + std::string(Foam::name(limitVal).c_str()) + " ]";
        }

        if (isBinding) line += "  <--";
        return line;
    };

    Info<< "  Time step control" << nl
        << "  " << std::string(52, '-').c_str() << nl
        << "  current deltaT:           " << actualDeltaT << nl;

    if (limitDielectricRelaxationRatio_ || printDielectricRelaxationRatio_)
    {
        const scalar val = (actualDeltaT * maxSigma) / eps0;
        const bool binding =
            limitDielectricRelaxationRatio_
        && (actualDeltaT * maxSigma) / eps0 >= maxDielectricRelaxationRatio_ * 0.99;

        Info<< "  " << fmtLine
            (
                "Diel. relax. ratio",
                val,
                limitDielectricRelaxationRatio_,
                maxDielectricRelaxationRatio_,
                binding
            ).c_str() << nl;
    }

    if (limitSpeciesCo_ || printSpeciesCo_)
    {
        const scalar convVal = maxConvFluxRate * actualDeltaT;
        const scalar diffVal = maxDiffFluxRate * actualDeltaT;
        const bool convBinding =
            limitSpeciesCo_ && convVal >= maxSpeciesConvectiveCo_ * 0.99;
        const bool diffBinding =
            limitSpeciesCo_ && diffVal >= maxSpeciesDiffusiveCo_ * 0.99;

        Info<< "  " << fmtLine
            (
                "Co_conv (" + courantSpeciesName_ + ")",
                convVal,
                limitSpeciesCo_,
                maxSpeciesConvectiveCo_,
                convBinding
            ).c_str() << nl
            << "  " << fmtLine
            (
                "Co_diff (" + courantSpeciesName_ + ")",
                diffVal,
                limitSpeciesCo_,
                maxSpeciesDiffusiveCo_,
                diffBinding
            ).c_str() << nl;
    }

    if (limitChemistryCo_ || printChemistryCo_)
    {
        const scalar val = maxKeff * actualDeltaT;
        const bool binding = limitChemistryCo_ && val >= maxChemistryCo_ * 0.99;

        Info<< "  " << fmtLine
            (
                "Co_chem",
                val,
                limitChemistryCo_,
                maxChemistryCo_,
                binding
            ).c_str() << nl;
    }

    if (limitVoltageRiseRate_ || printVoltageRiseRate_)
    {
        const scalar dVPerStep = voltageRiseRate * actualDeltaT;
        const bool binding =
            limitVoltageRiseRate_
        && dVPerStep >= maxVoltageRiseRate_ * 0.99;

        Info<< "  " << fmtLine
            (
                "dV/step [V]",
                dVPerStep,
                limitVoltageRiseRate_,
                maxVoltageRiseRate_,
                binding
            ).c_str() << nl;
    }

    Info<< "  " << std::string(52, '-').c_str() << nl << endl;
}

void plasmaTimeControl::setInitialDeltaT(const plasmaTransport& transport)
{
    if (!adjustTimeStep_) return;

    scalar newDeltaT = maxDeltaT_;
    const scalar currentDeltaT = runTime_.deltaTValue();
    const scalar eps0 = constant::plasma::epsilon0.value();

    // Dielectric relaxation (tau = epsilon / sigma)
    if (limitDielectricRelaxationRatio_)
    {
        tmp<volScalarField> tSigma = transport.electricalConductivity();
        const volScalarField& sigma = tSigma();
        const scalar maxSigma = gMax(mag(sigma)().primitiveField());
        newDeltaT = min
        (
            newDeltaT,
            (maxDielectricRelaxationRatio_ * eps0) / (maxSigma + VSMALL)
        );
        tSigma.clear();
    }

    // Species Courant Limit (convective and diffusive)
    if (limitSpeciesCo_)
    {
        const label sIdx = transport.species().speciesID(courantSpeciesName_);
        const volScalarField& n = transport.species().numberDensity(sIdx);
        const scalarField& nField = n.primitiveField();
        const scalarField& V = mesh_.V().field();

        // Convective
        {
            const surfaceScalarField& convFlux =
                transport.convectiveFlux(sIdx);
            const scalarField sumConv
            (
                fvc::surfaceSum(mag(convFlux))().primitiveField()
            );
            const scalar maxRate =
                gMax(0.5 * sumConv / (nField * V + VSMALL));
            newDeltaT = min
            (
                newDeltaT,
                maxSpeciesConvectiveCo_ / (maxRate + VSMALL)
            );
        }

        // Diffusive
        {
            const surfaceScalarField& diffFlux =
                transport.diffusiveFlux(sIdx);
            const scalarField sumDiff
            (
                fvc::surfaceSum(mag(diffFlux))().primitiveField()
            );
            const scalar maxRate =
                gMax(0.5 * sumDiff / (nField * V + VSMALL));
            newDeltaT = min
            (
                newDeltaT,
                maxSpeciesDiffusiveCo_ / (maxRate + VSMALL)
            );
        }
    }

    // Chemistry Courant Limit
    if (limitChemistryCo_)
    {
        // const volScalarField& keff = transport.k_eff();
        // const scalar maxKeff = gMax(mag(keff)().primitiveField());
        // newDeltaT = min
        // (
        //     newDeltaT,
        //     maxChemistryCo_ / (maxKeff + VSMALL)
        // );
    }

    // Voltage rise rate
    if (limitVoltageRiseRate_ || printVoltageRiseRate_)
    {
        prevPatchVoltage_ = patchVoltageAvg(transport);
    }

    // Reduction and setting
    reduce(newDeltaT, minOp<scalar>());
    if (newDeltaT > currentDeltaT)
    {
        newDeltaT = min(newDeltaT, currentDeltaT * 1.2);
    }
    runTime_.setDeltaT(newDeltaT);   
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
