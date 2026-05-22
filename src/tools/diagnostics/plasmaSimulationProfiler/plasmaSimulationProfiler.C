/*---------------------------------------------------------------------------*\
  File: plasmaSimulationProfiler.C
  Part of: SoPLASMA
  Developed using the OpenFOAM framework and linked against OpenFOAM libraries.

  Description:
    Implementation of Foam::plasmaSimulationProfiler.

  Copyright (C) 2026 Rention Pasolari
  License: GNU General Public License v3 or later
      See: <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "plasmaSimulationProfiler.H"
#include "IOstreams.H"
#include <iomanip>
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Initialize static members
    std::map<plasmaSimulationProfiler::ProfKey, double> plasmaSimulationProfiler::cumulativeTimes_;
    std::map<plasmaSimulationProfiler::ProfKey, double> plasmaSimulationProfiler::startTimes_;
    cpuTime plasmaSimulationProfiler::timer_;

    void plasmaSimulationProfiler::start(const std::string& main, const std::string& sub)
    {
        startTimes_[{main, sub}] = timer_.elapsedCpuTime();
    }

    void plasmaSimulationProfiler::stop(const std::string& main, const std::string& sub)
    {
        ProfKey key{main, sub};
        if (startTimes_.count(key))
        {
            cumulativeTimes_[key] += (timer_.elapsedCpuTime()
                                   - startTimes_[key]);
        }
    }

    void Foam::plasmaSimulationProfiler::report()
    {
        if (Pstream::master())
        {
            Info<< nl 
                << "===================== PROFILING REPORT ===================="
                << endl
                << "Section / Subsection           Max [s]   Avg [s]   Imbal %"
                << endl
                << "-----------------------------------------------------------"
                << endl;
        }

        for (auto const& [key, time] : cumulativeTimes_)
        {
            double maxT = time;
            double avgT = time;

            reduce(maxT, maxOp<double>());
            reduce(avgT, sumOp<double>());
            avgT /= Pstream::nProcs();

            double imbal = (maxT > 1e-6) ? (1.0 - (avgT / maxT)) * 100.0 : 0.0;

            if (Pstream::master())
            {
                std::string label = key.first;
                if (!key.second.empty()) label += " -> " + key.second;

                const int labelLimit = 50; // Increased from 30 to 50

                if (label.length() > labelLimit) 
                {
                    label = label.substr(0, labelLimit - 3) + "...";
                }
                else 
                {
                    label.append(labelLimit - label.length(), ' ');
                }

                // Use an ostringstream to build the formatted line
                std::ostringstream os;
                os << "\t" << label << " " 
                << std::fixed << std::setprecision(3) 
                << std::setw(9) << maxT << " "
                << std::setw(9) << avgT << " " 
                << std::setw(7) << imbal << "%";

                Info << os.str().c_str() << endl;
            }
        }

        if (Pstream::master())
        {
            Info<< "===========================================================" 
                << nl << endl;
        }
    }

    void plasmaSimulationProfiler::reset()
    {
        cumulativeTimes_.clear();
        startTimes_.clear();
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

 // End namespace Foam

// ************************************************************************* //
