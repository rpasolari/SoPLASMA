/*---------------------------------------------------------------------------*\
  File: plasmaTransportI.H
  Part of: foamPlasmaToolkit
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Inline member function definitions for the plasmaTransport class.

  Copyright (C) 2025 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef Foam_plasmaTransportI_H
#define Foam_plasmaTransportI_H

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Mesh and species metadata access
inline const Foam::fvMesh& 
Foam::plasmaTransport::mesh() const
{
    return mesh_;
}

inline const Foam::plasmaSpecies& 
Foam::plasmaTransport::species() const
{
    return species_;
}


//- Model access
inline Foam::plasmaTransportModel& 
Foam::plasmaTransport::model(const label i)
{
    return transportModels_[i];
}

inline const Foam::plasmaTransportModel&
Foam::plasmaTransport::model(const label i) const
{
    return transportModels_[i];
}

//- Total number flux access
inline const Foam::surfaceScalarField&
Foam::plasmaTransport::particleFlux(Foam::label i) const
{
    return particleFlux_[i];
}

inline Foam::surfaceScalarField&
Foam::plasmaTransport::particleFlux(Foam::label i)
{
    return particleFlux_[i];
}

inline const Foam::PtrList<Foam::surfaceScalarField>&
Foam::plasmaTransport::particleFlux() const
{
    return particleFlux_;
}

inline Foam::PtrList<Foam::surfaceScalarField>&
Foam::plasmaTransport::particleFlux()
{
    return particleFlux_;
}

//- Return the total volumetrix flux (m3/s) for a species on the fly
inline Foam::tmp<Foam::surfaceScalarField> 
Foam::plasmaTransport::phi(const label i) const
{
    return transportModels_[i].phi(); 
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
