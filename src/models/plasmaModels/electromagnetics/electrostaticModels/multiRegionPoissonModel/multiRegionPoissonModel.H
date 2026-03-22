#ifndef multiRegionPoissonModel_H
#define multiRegionPoissonModel_H

#include "electromagneticsModel.H"
#include "fvMatrix.H"

namespace Foam
{

class multiRegionPoissonModel
:
    public electromagneticsModel
{
    //--- Fields ---
    volScalarField ePotential_;
    volVectorField E_;
    dimensionedScalar epsilon_;

    PtrList<volScalarField> ePotentialDielectricList_;
    PtrList<volVectorField> EDielectricList_;
    PtrList<dimensionedScalar> epsilonDielectricList_;

    //--- Persistent Coupling State (Your Logic) ---
    bool coupled_;
    autoPtr<fvMatrix<scalar>> fvMatrixAssemblyPtr_;

public:

    TypeName("multiRegionPoissonModel");

    // FULL SIGNATURE
    multiRegionPoissonModel
    (
        const dynamicFvMesh& mesh,
        const UPtrList<dynamicFvMesh>& dielectricMeshes,
        const plasmaSpecies& species,
        const plasmaTransport& transport,
        const dictionary& dict
    );

    virtual ~multiRegionPoissonModel() {}

    virtual void solve() override;

    virtual const volVectorField& E() const override
    {
        return E_;
    }
};

} // End namespace Foam

#endif
