/*---------------------------------------------------------------------------*\
  File: plasmaTimeControl.C
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaTimeControl.

  Copyright (C) 2025 Rention Pasolari
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
    maxSpeciesCo_(1.0),
    speciesName_("e"),
    limitChemistryCo_(false),
    printChemistryCo_(false),
    maxChemistryCo_(1.0)
{
    read();
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

        limitDielectricRelaxationRatio_ =
         dict_.lookupOrDefault<Switch>("limitDielectricRelaxationRatio", false);

        printDielectricRelaxationRatio_ =
         dict_.lookupOrDefault<Switch>("printDielectricRelaxationRatio", false);

        maxDielectricRelaxationRatio_ =
            dict_.lookupOrDefault<scalar>("maxDielectricRelaxationRatio", 1.0);

        limitSpeciesCo_ = 
            dict_.lookupOrDefault<Switch>("limitSpeciesCo", false);

        printSpeciesCo_ = 
            dict_.lookupOrDefault<Switch>("printSpeciesCo", false);

        if (limitSpeciesCo_)
                {
                    maxSpeciesCo_ = 
                        dict_.lookupOrDefault<scalar>("maxSpeciesCo", 1.0);

                    speciesName_ = 
                        dict_.lookupOrDefault<word>("speciesName", "e");
                }

        limitChemistryCo_ = 
            dict_.lookupOrDefault<Switch>("limitChemistryCo", false);

        printChemistryCo_ = 
            dict_.lookupOrDefault<Switch>("printChemistryCo", false);

        maxChemistryCo_ = 
            dict_.lookupOrDefault<scalar>("maxChemistryCo", 1.0);
    }
}

void plasmaTimeControl::adjustDeltaT(const plasmaTransport& transport)
{
    const scalar currentDeltaT = runTime_.deltaTValue();
    const scalar eps0 = constant::plasma::epsilon0.value();
    scalar newDeltaT = maxDeltaT_;

    scalar maxSigma = 0.0;
    scalar maxFluxRate = 0.0;
    scalar meanFluxRate = 0.0;
    scalar maxKeff = 0.0;

    // Dielectric relaxation (tau = epsilon / sigma)
    if (limitDielectricRelaxationRatio_ || printDielectricRelaxationRatio_)
    {
        tmp<volScalarField> tSigma = transport.electricalConductivity();
        const volScalarField& sigma = tSigma();

        maxSigma = gMax(mag(sigma)().primitiveField());
        scalar dielectricLimit = 
            (maxDielectricRelaxationRatio_ * eps0) / (maxSigma + VSMALL);

        newDeltaT = min(newDeltaT, dielectricLimit);
        tSigma.clear();
    }

    // Species Courant Limit
    if (limitSpeciesCo_ || printSpeciesCo_)
    {
        label speciesLabel = transport.species().speciesID(speciesName_);
        const surfaceScalarField& phi = transport.phi(speciesLabel);

        scalarField sumPhi
        (
            fvc::surfaceSum(mag(phi))().primitiveField()
        );

        maxFluxRate = 0.5 * gMax(sumPhi / mesh_.V().field());
        meanFluxRate = 0.5 * (gSum(sumPhi) / gSum(mesh_.V().field()));

        scalar courantLimit = maxSpeciesCo_ / (maxFluxRate + VSMALL);

        newDeltaT = min(newDeltaT, courantLimit);
    }

    if (limitChemistryCo_ || printChemistryCo_)
    {
        // Access the k_eff field from transport
        const volScalarField& keff = transport.k_eff();
        
        maxKeff = gMax(mag(keff)().primitiveField());
        
        // dt <= maxChemistryCo / k_eff
        scalar chemistryLimit = maxChemistryCo_ / (maxKeff + VSMALL);

        newDeltaT = min(newDeltaT, chemistryLimit);
    }

    // Reduction and setting
    if (adjustTimeStep_)
    {
        // Global reduction for the actual time step decision
        reduce(newDeltaT, minOp<scalar>());

        if (newDeltaT > currentDeltaT)
        {
            newDeltaT = min(newDeltaT, currentDeltaT * 1.2);
        }
        runTime_.setDeltaT(newDeltaT);
    }

    // Report
    Info << "Time step monitoring:" << endl;
    Info << "  current deltaT   = " << currentDeltaT << endl;

    if (limitDielectricRelaxationRatio_ || printDielectricRelaxationRatio_)
    {
        scalar currentRatio = (currentDeltaT * maxSigma) / eps0;
        Info << "  deltaT/RelaxTime = " << currentRatio << endl;
    }

    if (limitSpeciesCo_ || printSpeciesCo_)
    {
        scalar currentCo = maxFluxRate * currentDeltaT;
        Info << "  Max Courant ("<< speciesName_ <<") = " << currentCo << endl;
    }

    if (limitChemistryCo_ || printChemistryCo_)
    {
        scalar currentChemCo = maxKeff * currentDeltaT;
        Info << "  Max Chem Courant = " << currentChemCo << endl;
    }
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

        scalar maxSigma = gMax(mag(sigma)().primitiveField());
        scalar dielectricLimit = 
            (maxDielectricRelaxationRatio_ * eps0) / (maxSigma + VSMALL);

        newDeltaT = min(newDeltaT, dielectricLimit);
        tSigma.clear();
    }

    // Species Courant Limit
    if (limitSpeciesCo_)
    {
        label speciesLabel = transport.species().speciesID(speciesName_);
        const surfaceScalarField& phi = transport.phi(speciesLabel);

        scalarField sumPhi
        (
            fvc::surfaceSum(mag(phi))().primitiveField()
        );

        scalar maxFluxRate = 0.5 * gMax(sumPhi / mesh_.V().field());
        scalar courantLimit = maxSpeciesCo_ / (maxFluxRate + VSMALL);

        newDeltaT = min(newDeltaT, courantLimit);
    }

    // Chemistry Courant Limit
    if (limitChemistryCo_)
    {
        const volScalarField& keff = transport.k_eff();
        
        scalar maxKeff = gMax(mag(keff)().primitiveField());
        scalar chemistryLimit = maxChemistryCo_ / (maxKeff + VSMALL);

        newDeltaT = min(newDeltaT, chemistryLimit);
    }

    // Reduction and setting
    reduce(newDeltaT, minOp<scalar>());
    if (newDeltaT < currentDeltaT)
    {
        runTime_.setDeltaT(newDeltaT);
    }    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
