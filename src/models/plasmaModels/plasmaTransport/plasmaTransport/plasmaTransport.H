/*---------------------------------------------------------------------------*\
License
    This file is part of the foamPlasmaToolkit.

    The foamPlasmaToolkit is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2025
        Rention Pasolari

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::plasmaTransport

Description
    High-level manager for plasma transport models.

    This class owns and manages one transport model per plasma species.
    It does not implement any transport physics itself. Instead, it serves
    as a controller that:

        - Retrieves species-level transport configuration from plasmaSpecies
        - Constructs the appropriate plasmaTransportModel for each species
        - Stores all per-species transport models in a PtrList
        - Delegates time-dependent updates to the individual models

    The plasmaTransport object is typically created once in the solver and
    used during the solution loop to update the transport behaviour of all
    plasma species.

SourceFiles
    plasmaTransportI.H
    plasmaTransport.C
    
\*---------------------------------------------------------------------------*/

#ifndef Foam_plasmaTransport_H
#define Foam_plasmaTransport_H

#include "fvCFD.H"
#include "fvMesh.H"
#include "PtrList.H"

#include "plasmaSpecies.H"
#include "plasmaTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class plasmaTransport Declaration
\*---------------------------------------------------------------------------*/

class plasmaTransport
:
    public regIOobject
{
private:

    //- Reference to mesh
    const fvMesh& mesh_;

    //- Reference to species database
    plasmaSpecies& species_;

    //- Reference to the electric field
    const volVectorField& E_;

    //- Reference to the electric field flux
    surfaceScalarField phiE_;

    //- Reference to the electric potential field
    const volScalarField& ePotential_;

    //- Total particle surface Flux [1/s]
    PtrList<surfaceScalarField> particleFlux_;

    //- Effective reaction rate field [1/s]
    volScalarField k_eff_;

    //- Per-species transport models
    PtrList<plasmaTransportModel> transportModels_;

    //- Disallow copying
    plasmaTransport(const plasmaTransport&) = delete;
    void operator=(const plasmaTransport&) = delete;

    //- Construct per-species models
    void constructModels();

public:

    //- Runtime type information
    TypeName("plasmaTransport");

    // Constructors

        //- Construct from plasmaSpecies and mesh
        plasmaTransport
        (
            plasmaSpecies& species,
            const fvMesh& mesh,
            const volVectorField& E,
            const volScalarField& ePotential
        );

    //- Destructor
    virtual ~plasmaTransport() {}

    // Member functions

        //- Correct all models
        void correct(const bool firstIter, const bool finalIter);

        //- Return total electrical conductivity [S/m]
        tmp<volScalarField> electricalConductivity() const;

        //- Return total diffusive charge source [C/m^3/s]
        tmp<volScalarField> diffusiveChargeSource() const;

        virtual bool writeData(Ostream& os) const;

    // Access

        //- Mesh and species metadata access
        inline const fvMesh& mesh() const;
        inline const plasmaSpecies& species() const;

        //- Model access
        inline plasmaTransportModel& model(const label i);
        inline const plasmaTransportModel& model(const label i) const;

        //- Total number flux access
        inline const surfaceScalarField& particleFlux(label i) const;
        inline surfaceScalarField& particleFlux(label i);
        inline const PtrList<surfaceScalarField>& particleFlux() const;
        inline PtrList<surfaceScalarField>& particleFlux();

        //- Return the effective reaction rate field [1/s]
        inline const volScalarField& k_eff() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "plasmaTransportI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
