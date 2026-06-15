/*---------------------------------------------------------------------------*\
License
    This file is part of the SoPLASMA.

    The SoPLASMA is not part of OpenFOAM but is developed using the
    OpenFOAM framework and linked against OpenFOAM libraries.

    Copyright (C) 2026
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
    Manager class for per-species plasma transport models.

    Owns one plasmaTransportModel per active species. Delegates all
    transport physics to the individual models and aggregates their
    outputs for use by the electromagnetics solver.

    Provides:
        - Per-species particle flux (for number density equations)
        - Total electrical conductivity (for semi-implicit Poisson)
        - Total diffusive charge source (for semi-implicit Poisson)
        - Wall flux list (for surface charge update via plasmaSpecies)

    Has no knowledge of reactions or energy models.

SourceFiles
    plasmaTransportI.H
    plasmaTransport.C
    
\*---------------------------------------------------------------------------*/

#ifndef Foam_plasmaTransport_H
#define Foam_plasmaTransport_H

#include "fvCFD.H"
#include "fvMesh.H"
#include "PtrList.H"
#include "regIOobject.H"

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

    //- Disallow copy
    plasmaTransport(const plasmaTransport&) = delete;
    void operator=(const plasmaTransport&) = delete;

    //- Construct per-species transport models
    void constructTransportModels();

protected:

    //- Reference to mesh
    const fvMesh& mesh_;

    //- Reference to species database
    plasmaSpecies& species_;

    //- Per-species transport models
    PtrList<plasmaTransportModel> transportModels_;

    //- Convective particle flux [1/s]
    PtrList<surfaceScalarField> convectiveFlux_;

    //- Diffusive particle flux [1/s]
    PtrList<surfaceScalarField> diffusiveFlux_;

    //- Total particle surface flux [1/s]
    PtrList<surfaceScalarField> particleFlux_;

    volScalarField k_eff_;
    volScalarField alpha_;
     volScalarField alphaDx_;

public:

    //- Runtime type information
    TypeName("plasmaTransport");

    // Constructors

        //- Construct from mesh and species
        plasmaTransport
        (
            const fvMesh& mesh,
            plasmaSpecies& species
        );

    //- Destructor
    virtual ~plasmaTransport() = default;

    // Member functions

        //- Correct all transport models
        void correctTransportModels();

        //- Solve all species transport equations
        void solve();

        //- Update surface charge
        void updateSurfaceCharge();

        //- Return the LHS matrix for a specific species
        tmp<fvScalarMatrix> transportMatrix(const label i);

        //- Return total electrical conductivity [S/m]
        tmp<volScalarField> electricalConductivity() const;

        //- Return total diffusive charge source [C/m^3/s]
        tmp<volScalarField> diffusiveChargeSource() const;

        virtual bool writeData(Ostream& os) const;

    // Access

        //- Metadata access
        inline const fvMesh& mesh() const;
        inline plasmaSpecies& species();
        inline const plasmaSpecies& species() const;

        //- Access the EM model through species
        inline const electromagneticsModel& em() const;

        //- Model access
        inline plasmaTransportModel& model(const label i);
        inline const plasmaTransportModel& model(const label i) const;

        //- Access fluxes
        inline surfaceScalarField& particleFlux(const label i);
        inline const surfaceScalarField& particleFlux(const label i) const;
        inline PtrList<surfaceScalarField>& particleFlux();
        inline const PtrList<surfaceScalarField>& particleFlux() const;
        inline surfaceScalarField& convectiveFlux(const label i);
        inline const surfaceScalarField& convectiveFlux(const label i) const;
        inline PtrList<surfaceScalarField>& convectiveFlux();
        inline const PtrList<surfaceScalarField>& convectiveFlux() const;
        inline surfaceScalarField& diffusiveFlux(const label i);
        inline const surfaceScalarField& diffusiveFlux(const label i) const;
        inline PtrList<surfaceScalarField>& diffusiveFlux();
        inline const PtrList<surfaceScalarField>& diffusiveFlux() const;

        inline const volScalarField& k_eff() const { return k_eff_; };
        const volScalarField& alpha() const noexcept { return alpha_; }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "plasmaTransportI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
