#ifndef electromagneticsModel_H
#define electromagneticsModel_H

#include "dynamicFvMesh.H"
#include "UPtrList.H"
#include "runTimeSelectionTables.H"
#include "plasmaSpecies.H"   // Based on your code: plasmaSpecies
#include "plasmaTransport.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class electromagneticsModel Declaration
\*---------------------------------------------------------------------------*/

class electromagneticsModel
{
protected:

    //- Reference to the primary gas/plasma dynamic mesh
    const dynamicFvMesh& mesh_;

    //- Reference to the list of dielectric dynamic meshes (non-owning)
    const UPtrList<dynamicFvMesh>& dielectricMeshes_;

    //- Reference to plasma composition/species
    const plasmaSpecies& species_;

    //- Reference to plasma transport models
    const plasmaTransport& transport_;

    //- Reference to the electromagnetics sub-dictionary
    const dictionary& dict_;


public:

    //- Runtime type information
    TypeName("electromagneticsModel");

    //- Declare the Run-Time Selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        electromagneticsModel,
        dictionary,
        (
            const dynamicFvMesh& mesh,
            const UPtrList<dynamicFvMesh>& dielectricMeshes,
            const plasmaSpecies& species,
            const plasmaTransport& transport,
            const dictionary& dict
        ),
        (mesh, dielectricMeshes, species, transport, dict)
    );


    //--- Constructors ---

    electromagneticsModel
    (
        const dynamicFvMesh& mesh,
        const UPtrList<dynamicFvMesh>& dielectricMeshes,
        const plasmaSpecies& species,
        const plasmaTransport& transport,
        const dictionary& dict
    );


    //- Selector: The "Factory" function that reads the dictionary
    static autoPtr<electromagneticsModel> New
    (
        const dynamicFvMesh& mesh,
        const UPtrList<dynamicFvMesh>& dielectricMeshes,
        const plasmaSpecies& species,
        const plasmaTransport& transport,
        const dictionary& dict
    );


    //--- Destructor ---

    virtual ~electromagneticsModel()
    {}


    //--- Member Functions ---

    //- The main solve command (must be implemented by children)
    virtual void solve() = 0;

    // Inside electromagneticsModel.H
virtual const volVectorField& E() const = 0;

    //- Access to the main mesh
    const dynamicFvMesh& mesh() const
    {
        return mesh_;
    }
};

} // End namespace Foam

#endif
