#ifndef singleRegionPoissonModel_H
#define singleRegionPoissonModel_H

#include "electromagneticsModel.H"

namespace Foam
{

class singleRegionPoissonModel
:
    public electromagneticsModel
{
    //- Fields owned by this class
    volScalarField ePotential_;
    volVectorField E_;

    //- Permittivity
    dimensionedScalar epsilon_;

public:

    TypeName("singleRegionPoisson");

    singleRegionPoissonModel
    (
        const dynamicFvMesh& mesh,
        const UPtrList<dynamicFvMesh>& dielectricMeshes,
        const plasmaSpecies& species,
        const plasmaTransport& transport,
        const dictionary& dict
    );

    virtual ~singleRegionPoissonModel() {}

    virtual void solve() override;

    //- Getter for the Electric Field
    virtual const volVectorField& E() const override;
};

} // End namespace Foam

#endif
