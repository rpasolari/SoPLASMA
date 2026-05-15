/*---------------------------------------------------------------------------*\
  File: createDielectricFields.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Derived from createSolidFields.H in OpenFOAM v2412.
      Defines and initializes the dielectric-region fields for each region.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

PtrList<volScalarField> ePotentialDielectric(dielectricRegions.size());
PtrList<IOdictionary> electricPropertiesDictDielectric(dielectricRegions.size());
PtrList<volScalarField> surfChargeDielectric(dielectricRegions.size());
PtrList<dimensionedScalar> epsilonDielectric(dielectricRegions.size());

List<scalar> epsilonRDielectric(dielectricRegions.size());

List<scalar> ePotentialNonCoupledResidualDielectric
(
    dielectricRegions.size(), Foam::GREAT
);

forAll(dielectricRegions, i)
{
    Info << "*** Reading fields for region "
         << dielectricRegions[i].name() << nl << endl;

    Info << "   Adding to ePotentialDielectric\n" << endl;
    ePotentialDielectric.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "ePotential",
                runTime.timeName(),
                dielectricRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            dielectricRegions[i]
        )
    );

    Info << "   Adding to surfChargeDielectric\n" << endl;
    surfChargeDielectric.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "surfCharge",
                runTime.timeName(),
                dielectricRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            dielectricRegions[i],
            dimensionedScalar("zero", dimensionSet(0, -2, 1, 0, 0, 1, 0), 0.0)
        )
    );

    electricPropertiesDictDielectric.set
    (
        i,
        new IOdictionary
        (
            IOobject
            (
                "electricProperties",
                runTime.constant(),
                dielectricRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    epsilonRDielectric[i] =
        electricPropertiesDictDielectric[i].lookupOrDefault<scalar>(
            "dielectricConstant", 1.0
        );

    Info << "   Adding material's permittivity\n" << endl;
    epsilonDielectric.set
    (
        i,
        new dimensionedScalar
        (
            "epsilonDielectric",
            epsilonRDielectric[i] * constant::plasma::epsilon0
        )
    );
}
