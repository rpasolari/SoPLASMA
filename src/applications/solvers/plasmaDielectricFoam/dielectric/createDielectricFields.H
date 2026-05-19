/*---------------------------------------------------------------------------*\
  File: createDielectricFields.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Derived from createSolidFields.H in OpenFOAM v2412.
      Defines and initializes the dielectric-region fields for each region.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

PtrList<volScalarField> ePotentialDielectricList(dielectricRegions.size());
PtrList<volVectorField> EDielectricList(dielectricRegions.size());
PtrList<IOdictionary>
    electricPropertiesDictDielectricList(dielectricRegions.size());
PtrList<volScalarField> surfChargeDielectricList(dielectricRegions.size());
PtrList<dimensionedScalar> epsilonDielectricList(dielectricRegions.size());

List<scalar> epsilonRDielectricList(dielectricRegions.size());

List<scalar> ePotentialNonCoupledResidualDielectricList
(
    dielectricRegions.size(), Foam::GREAT
);

forAll(dielectricRegions, i)
{
    Info << "*** Reading fields for region "
         << dielectricRegions[i].name() << nl << endl;

    Info << "   Adding to ePotentialDielectric\n" << endl;
    ePotentialDielectricList.set
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

    Info << "   Adding to EDielectric\n" << endl;
    EDielectricList.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "E",
                runTime.timeName(),
                dielectricRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            dielectricRegions[i],
            dimensionedVector
            (
                "zero",
                dimensionSet(1, 1, -3, 0, 0, -1, 0),
                vector::zero
            )
        )
    );

    Info << "   Adding to surfChargeDielectric\n" << endl;
    surfChargeDielectricList.set
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

    electricPropertiesDictDielectricList.set
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

    if (!electricPropertiesDictDielectricList[i].found("dielectricConstant"))
    {
        FatalErrorInFunction
            << "Missing 'dielectricConstant' in electricProperties for region "
            << dielectricRegions[i].name()
            << "\n    Please specify the relative permittivity."
            << exit(FatalError);
    }

    epsilonRDielectricList[i] =
        electricPropertiesDictDielectricList[i].get<scalar>
        ("dielectricConstant");

    Info << "   Adding material's permittivity\n" << endl;
    epsilonDielectricList.set
    (
        i,
        new dimensionedScalar
        (
            "epsilonDielectric",
            epsilonRDielectricList[i] * constant::plasma::epsilon0
        )
    );
}
