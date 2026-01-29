/*---------------------------------------------------------------------------*\
License
    This file is part of the foamPlasmaToolkit.

    The foamPlasmaToolkit is not part of OpenFOAM but is developed using the
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
    Foam::plasmaProfiler

Description

SourceFiles
    plasmaProfiler.C

\*---------------------------------------------------------------------------*/

#ifndef plasmaProfiler_H
#define plasmaProfiler_H

#include "cpuTime.H"
#include "Pstream.H"
#include <map>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class plasmaProfiler Declaration
\*---------------------------------------------------------------------------*/

class plasmaProfiler
{
private:

    // Private data

        //- Unique ID: {MainCategory, SubCategory}
        typedef std::pair<std::string, std::string> ProfKey;

        //- Total seconds accumulated for each ID
        static std::map<ProfKey, double> cumulativeTimes_;

        //- Starting timestamp for active measurements
        static std::map<ProfKey, double> startTimes_;

        //- The master stopwatch for the library
        static cpuTime timer_;
    
public:

    //- Start timing a section. Sub-tag is optional
    static void start(const std::string& main, const std::string& sub = "");

    //- Stop timing and save the duration
    static void stop(const std::string& main, const std::string& sub = "");

    //- Synchronize parallel data and print the report
    static void report();

    //- Clear all stored timing data
    static void reset();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
