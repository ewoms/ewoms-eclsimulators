/*

  This file is part of the eWoms project.

  eWoms is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  eWoms is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with eWoms.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "config.h"

#include <ewoms/eclsimulators/timestepping/gatherconvergencereport.hh>

#if HAVE_MPI

#include <mpi.h>

namespace
{

    using Ewoms::ConvergenceReport;

    void packReservoirFailure(const ConvergenceReport::ReservoirFailure& f,
                              std::vector<char>& buf,
                              int& offset)
    {
        int type = static_cast<int>(f.type());
        int severity = static_cast<int>(f.severity());
        int phase = f.phase();
        MPI_Pack(&type, 1, MPI_INT, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
        MPI_Pack(&severity, 1, MPI_INT, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
        MPI_Pack(&phase, 1, MPI_INT, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
    }

    void packWellFailure(const ConvergenceReport::WellFailure& f,
                         std::vector<char>& buf,
                         int& offset)
    {
        int type = static_cast<int>(f.type());
        int severity = static_cast<int>(f.severity());
        int phase = f.phase();
        MPI_Pack(&type, 1, MPI_INT, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
        MPI_Pack(&severity, 1, MPI_INT, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
        MPI_Pack(&phase, 1, MPI_INT, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
        int name_length = f.wellName().size() + 1; // Adding 1 for the null terminator.
        MPI_Pack(&name_length, 1, MPI_INT, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
        MPI_Pack(const_cast<char*>(f.wellName().c_str()), name_length, MPI_CHAR, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
    }

    void packConvergenceReport(const ConvergenceReport& local_report,
                               std::vector<char>& buf,
                               int& offset)
    {
        // Pack the data.
        // Status will not be packed, it is possible to deduce from the other data.
        // Reservoir failures.
        const auto rf = local_report.reservoirFailures();
        int num_rf = rf.size();
        MPI_Pack(&num_rf, 1, MPI_INT, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
        for (const auto& f : rf) {
            packReservoirFailure(f, buf, offset);
        }
        // Well failures.
        const auto wf = local_report.wellFailures();
        int num_wf = wf.size();
        MPI_Pack(&num_wf, 1, MPI_INT, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
        for (const auto& f : wf) {
            packWellFailure(f, buf, offset);
        }
    }

    int messageSize(const ConvergenceReport& local_report)
    {
        int intPackSize = 0;
        MPI_Pack_size(1, MPI_INT, MPI_COMM_WORLD, &intPackSize);
        const int num_rf = local_report.reservoirFailures().size();
        const int num_wf = local_report.wellFailures().size();
        int wellnames_length = 0;
        for (const auto& f : local_report.wellFailures()) {
            wellnames_length += (f.wellName().size() + 1);
        }
        return (2 + 3*num_rf + 4*num_wf) * intPackSize + wellnames_length;
    }

    ConvergenceReport::ReservoirFailure unpackReservoirFailure(const std::vector<char>& recvBuffer, int& offset)
    {
        int type = -1;
        int severity = -1;
        int phase = -1;
        auto* data = const_cast<char*>(recvBuffer.data());
        MPI_Unpack(data, recvBuffer.size(), &offset, &type, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack(data, recvBuffer.size(), &offset, &severity, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack(data, recvBuffer.size(), &offset, &phase, 1, MPI_INT, MPI_COMM_WORLD);
        return ConvergenceReport::ReservoirFailure(static_cast<ConvergenceReport::ReservoirFailure::Type>(type),
                                                   static_cast<ConvergenceReport::Severity>(severity),
                                                   phase);
    }

    ConvergenceReport::WellFailure unpackWellFailure(const std::vector<char>& recvBuffer, int& offset)
    {
        int type = -1;
        int severity = -1;
        int phase = -1;
        auto* data = const_cast<char*>(recvBuffer.data());
        MPI_Unpack(data, recvBuffer.size(), &offset, &type, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack(data, recvBuffer.size(), &offset, &severity, 1, MPI_INT, MPI_COMM_WORLD);
        MPI_Unpack(data, recvBuffer.size(), &offset, &phase, 1, MPI_INT, MPI_COMM_WORLD);
        int name_length = -1;
        MPI_Unpack(data, recvBuffer.size(), &offset, &name_length, 1, MPI_INT, MPI_COMM_WORLD);
        std::vector<char> namechars(name_length);
        MPI_Unpack(data, recvBuffer.size(), &offset, namechars.data(), name_length, MPI_CHAR, MPI_COMM_WORLD);
        std::string name(namechars.data());
        return ConvergenceReport::WellFailure(static_cast<ConvergenceReport::WellFailure::Type>(type),
                                              static_cast<ConvergenceReport::Severity>(severity),
                                              phase,
                                              name);
    }

    ConvergenceReport unpackSingleConvergenceReport(const std::vector<char>& recvBuffer, int& offset)
    {
        ConvergenceReport cr;
        int num_rf = -1;
        auto* data = const_cast<char*>(recvBuffer.data());
        MPI_Unpack(data, recvBuffer.size(), &offset, &num_rf, 1, MPI_INT, MPI_COMM_WORLD);
        for (int rf = 0; rf < num_rf; ++rf) {
            ConvergenceReport::ReservoirFailure f = unpackReservoirFailure(recvBuffer, offset);
            cr.setReservoirFailed(f);
        }
        int num_wf = -1;
        MPI_Unpack(data, recvBuffer.size(), &offset, &num_wf, 1, MPI_INT, MPI_COMM_WORLD);
        for (int wf = 0; wf < num_wf; ++wf) {
            ConvergenceReport::WellFailure f = unpackWellFailure(recvBuffer, offset);
            cr.setWellFailed(f);
        }
        return cr;
    }

    ConvergenceReport unpackConvergenceReports(const std::vector<char>& recvBuffer,
                                               const std::vector<int>& displ)
    {
        ConvergenceReport cr;
        const int numProcesses = displ.size() - 1;
        for (int process = 0; process < numProcesses; ++process) {
            int offset = displ[process];
            cr += unpackSingleConvergenceReport(recvBuffer, offset);
            assert(offset == displ[process + 1]);
        }
        return cr;
    }

} // anonymous namespace

namespace Ewoms
{

    /// Create a global convergence report combining local
    /// (per-process) reports.
    ConvergenceReport gatherConvergenceReport(const ConvergenceReport& local_report)
    {
        // Pack local report.
        int messageSize = ::messageSize(local_report);
        std::vector<char> buffer(messageSize);
        int offset = 0;
        packConvergenceReport(local_report, buffer, offset);
        assert(offset == messageSize);

        // Get message sizes and create offset/displacement array for gathering.
        int numProcesses = -1;
        MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
        std::vector<int> messageSizes(numProcesses);
        MPI_Allgather(&messageSize, 1, MPI_INT, messageSizes.data(), 1, MPI_INT, MPI_COMM_WORLD);
        std::vector<int> displ(numProcesses + 1, 0);
        std::partial_sum(messageSizes.begin(), messageSizes.end(), displ.begin() + 1);

        // Gather.
        std::vector<char> recvBuffer(displ.back());
        MPI_Allgatherv(buffer.data(), buffer.size(), MPI_PACKED,
                       const_cast<char*>(recvBuffer.data()), messageSizes.data(),
                       displ.data(), MPI_PACKED,
                       MPI_COMM_WORLD);

        // Unpack.
        ConvergenceReport global_report = unpackConvergenceReports(recvBuffer, displ);
        return global_report;
    }

} // namespace Ewoms

#else // HAVE_MPI

namespace Ewoms
{
    ConvergenceReport gatherConvergenceReport(const ConvergenceReport& local_report)
    {
        return local_report;
    }
} // namespace Ewoms

#endif // HAVE_MPI
