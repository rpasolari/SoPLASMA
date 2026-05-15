/*---------------------------------------------------------------------------*\
  File: createFluidFields.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Derived from createFluidFields.H in OpenFOAM v2412.
      Defines and initializes the fluid-region fields for each region.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

PtrList<volScalarField> ePotentialFluid(fluidRegions.size());
PtrList<volScalarField> eChargeFluid(fluidRegions.size());
PtrList<volScalarField> surfChargeFluid(fluidRegions.size());
PtrList<IOdictionary> electricPropertiesDictFluid(fluidRegions.size());
PtrList<dimensionedScalar> epsilonFluid(fluidRegions.size());

List<scalar> epsilonRFluid(fluidRegions.size());

List<scalar> ePotentialNonCoupledResidualFluid
(
    fluidRegions.size(), Foam::GREAT
);

forAll(fluidRegions, i)
{
    Info << "*** Reading fields for region "
         << fluidRegions[i].name() << nl << endl;

    Info << "   Adding to ePotentialFluid\n" << endl;
    ePotentialFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "ePotential",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info << "   Adding to eChargeFluid\n" << endl;
    eChargeFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "eCharge",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info << "   Adding to surfChargeFluid\n" << endl;
    surfChargeFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "surfCharge",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar("zero", dimensionSet(0, -2, 1, 0, 0, 1, 0), 0.0)
        )
    );

    electricPropertiesDictFluid.set
    (
        i,
        new IOdictionary
        (
            IOobject
            (
                "electricProperties",
                runTime.constant(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    epsilonRFluid[i] = electricPropertiesDictFluid[i].lookupOrDefault<scalar>
    (
        "dielectricConstant", 1.0
    );

    Info << "   Adding material's permittivity\n" << endl;

    epsilonFluid.set
    (
        i,
        new dimensionedScalar
        (
            "epsilonFluid",
            epsilonRFluid[i] * constant::plasma::epsilon0
        )
    );
}


