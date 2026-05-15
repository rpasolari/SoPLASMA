/*---------------------------------------------------------------------------*\
  File: createFields.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Derived from createFields.H in OpenFOAM v2412.
      Defines and initializes electrostatic fields.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

Info << "*** Reading fields" << endl;

// Electrostatic potential
Info << "   Reading ePotential\n" << endl;
volScalarField ePotential
(
    IOobject
    (
        "ePotential",
        runTime.timeName(),
        mesh,                          
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Charge density
Info << "   Reading to eCharge\n" << endl;
volScalarField eCharge
(
    IOobject
    (
        "eCharge",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Read material electric properties
Info << "   Reading material's permittivity\n" << endl;
IOdictionary electricPropertiesDict
(
    IOobject
    (
        "electricProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

// Relative permittivity
scalar epsilonR =
    electricPropertiesDict.lookupOrDefault<scalar>("dielectricConstant", 1.0);

Info<< "   Material relative permittivity: " << epsilonR << nl;

// Absolute permittivity
dimensionedScalar epsilon
(
    "epsilon",
    epsilonR * constant::plasma::epsilon0
);

// ePotentialControls
dictionary epCtrl =
    mesh.solutionDict().subDict("ePotentialControls");

// non-orthogonal correctors
int nNonOrthoCorr =
    epCtrl.getOrDefault<int>("nNonOrthogonalCorrectors", 0);
