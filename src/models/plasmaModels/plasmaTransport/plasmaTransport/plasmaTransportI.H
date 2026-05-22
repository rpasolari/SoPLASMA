/*---------------------------------------------------------------------------*\
  File: plasmaTransportI.H
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
      Inline member function definitions for the plasmaTransport class.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#ifndef Foam_plasmaTransportI_H
#define Foam_plasmaTransportI_H

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Metadata access
inline const Foam::fvMesh&
Foam::plasmaTransport::mesh() const
{
    return mesh_;
}

inline Foam::plasmaSpecies&
Foam::plasmaTransport::species()
{
    return species_;
}

inline const Foam::plasmaSpecies&
Foam::plasmaTransport::species() const
{
    return species_;
}

inline const Foam::electromagneticsModel&
Foam::plasmaTransport::em() const
{
    return species_.em();
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

//- Fluxes access
inline Foam::surfaceScalarField&
Foam::plasmaTransport::particleFlux(const Foam::label i)
{
    return particleFlux_[i];
}

inline const Foam::surfaceScalarField&
Foam::plasmaTransport::particleFlux(const Foam::label i) const
{
    return particleFlux_[i];
}

inline Foam::PtrList<Foam::surfaceScalarField>&
Foam::plasmaTransport::particleFlux()
{
    return particleFlux_;
}

inline const Foam::PtrList<Foam::surfaceScalarField>&
Foam::plasmaTransport::particleFlux() const
{
    return particleFlux_;
}

inline Foam::surfaceScalarField&
Foam::plasmaTransport::convectiveFlux(const Foam::label i)
{
    return convectiveFlux_[i];
}

inline const Foam::surfaceScalarField&
Foam::plasmaTransport::convectiveFlux(const Foam::label i) const
{
    return convectiveFlux_[i];
}

inline Foam::PtrList<Foam::surfaceScalarField>&
Foam::plasmaTransport::convectiveFlux()
{
    return convectiveFlux_;
}

inline const Foam::PtrList<Foam::surfaceScalarField>&
Foam::plasmaTransport::convectiveFlux() const
{
    return convectiveFlux_;
}

inline Foam::surfaceScalarField&
Foam::plasmaTransport::diffusiveFlux(const Foam::label i)
{
    return diffusiveFlux_[i];
}

inline const Foam::surfaceScalarField&
Foam::plasmaTransport::diffusiveFlux(const Foam::label i) const
{
    return diffusiveFlux_[i];
}

inline Foam::PtrList<Foam::surfaceScalarField>&
Foam::plasmaTransport::diffusiveFlux()
{
    return diffusiveFlux_;
}

inline const Foam::PtrList<Foam::surfaceScalarField>&
Foam::plasmaTransport::diffusiveFlux() const
{
    return diffusiveFlux_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
