#include "electromagneticsModel.H"

namespace Foam
{
    // Define the Run-Time Selection table
    defineRunTimeSelectionTable(electromagneticsModel, dictionary);

    // Constructor implementation
    electromagneticsModel::electromagneticsModel
    (
        const dynamicFvMesh& mesh,
        const UPtrList<dynamicFvMesh>& dielectricMeshes,
        const plasmaSpecies& species,
        const plasmaTransport& transport,
        const dictionary& dict
    )
    :
        mesh_(mesh),
        dielectricMeshes_(dielectricMeshes),
        species_(species),
        transport_(transport),
        dict_(dict)
    {}


    // The "New" Selector (The Factory)
    autoPtr<electromagneticsModel> electromagneticsModel::New
    (
        const dynamicFvMesh& mesh,
        const UPtrList<dynamicFvMesh>& dielectricMeshes,
        const plasmaSpecies& species,
        const plasmaTransport& transport,
        const dictionary& dict
    )
    {
        // Look up the "model" keyword (e.g., singleRegionPoisson)
        const word modelType(dict.get<word>("model"));

        Info<< "Selecting electromagneticsModel " << modelType << endl;

        // Look for the model in the RTS table
        auto* ctorPtr = dictionaryConstructorTable(modelType);

        if (!ctorPtr)
        {
            FatalErrorInFunction
                << "Unknown electromagneticsModel type " << modelType
                << "\n\nValid models are: " 
                << dictionaryConstructorTablePtr_->sortedToc()
                << exit(FatalError);
        }

        // Create and return the chosen child class
        return autoPtr<electromagneticsModel>
        (
            ctorPtr(mesh, dielectricMeshes, species, transport, dict)
        );
    }
}
