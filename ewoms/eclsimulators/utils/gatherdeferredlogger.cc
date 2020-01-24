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

#include <ewoms/eclsimulators/utils/gatherdeferredlogger.hh>

#if HAVE_MPI

#include <cassert>
#include <cstdint>
#include <numeric>
#include <mpi.h>

namespace
{

    void packMessages(const std::vector<Ewoms::DeferredLogger::Message>& local_messages, std::vector<char>& buf, int& offset)
    {

        int messagesize = local_messages.size();
        MPI_Pack(&messagesize, 1, MPI_UNSIGNED, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);

        for (const auto lm : local_messages) {
            MPI_Pack(static_cast<void*>(const_cast<std::int64_t*>(&lm.flag)), 1, MPI_INT64_T, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
            int tagsize = lm.tag.size();
            MPI_Pack(&tagsize, 1, MPI_UNSIGNED, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
            if (tagsize>0) {
                MPI_Pack(const_cast<char*>(lm.tag.c_str()), lm.tag.size(), MPI_CHAR, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
            }
            int textsize = lm.text.size();
            MPI_Pack(&textsize, 1, MPI_UNSIGNED, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
            if (textsize>0) {
                MPI_Pack(const_cast<char*>(lm.text.c_str()), lm.text.size(), MPI_CHAR, buf.data(), buf.size(), &offset, MPI_COMM_WORLD);
            }
        }
    }

    Ewoms::DeferredLogger::Message unpackSingleMessage(const std::vector<char>& recvBuffer, int& offset)
    {
        int64_t flag;
        auto* data = const_cast<char*>(recvBuffer.data());
        MPI_Unpack(data, recvBuffer.size(), &offset, &flag, 1, MPI_INT64_T, MPI_COMM_WORLD);

        // unpack tag
        unsigned int tagsize;
        MPI_Unpack(data, recvBuffer.size(), &offset, &tagsize, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
        std::string tag;
        if (tagsize>0) {
            std::vector<char> tagchars(tagsize);
            MPI_Unpack(data, recvBuffer.size(), &offset, tagchars.data(), tagsize, MPI_CHAR, MPI_COMM_WORLD);
            tag = std::string(tagchars.data(), tagsize);
        }
        // unpack text
        unsigned int textsize;
        MPI_Unpack(data, recvBuffer.size(), &offset, &textsize, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
        std::string text;
        if (textsize>0) {
            std::vector<char> textchars(textsize);
            MPI_Unpack(data, recvBuffer.size(), &offset, textchars.data(), textsize, MPI_CHAR, MPI_COMM_WORLD);
            text = std::string (textchars.data(), textsize);
        }
        return Ewoms::DeferredLogger::Message({flag, tag, text});
    }

    std::vector<Ewoms::DeferredLogger::Message> unpackMessages(const std::vector<char>& recvBuffer, const std::vector<int>& displ)
    {
        std::vector<Ewoms::DeferredLogger::Message> messages;
        const int numProcesses = displ.size() - 1;
        auto* data = const_cast<char*>(recvBuffer.data());
        for (int process = 0; process < numProcesses; ++process) {
            int offset = displ[process];
            // unpack number of messages
            unsigned int messagesize;
            MPI_Unpack(data, recvBuffer.size(), &offset, &messagesize, 1, MPI_UNSIGNED, MPI_COMM_WORLD);
            for (unsigned int i=0; i<messagesize; i++) {
                messages.push_back(unpackSingleMessage(recvBuffer, offset));
            }
            assert(offset == displ[process + 1]);
        }
        return messages;
    }

} // anonymous namespace

namespace Ewoms
{

    /// combine (per-process) messages
    Ewoms::DeferredLogger gatherDeferredLogger(const Ewoms::DeferredLogger& local_deferredlogger)
    {

        int num_messages = local_deferredlogger.messages_.size();

        int int64MpiPackSize;
        MPI_Pack_size(1, MPI_INT64_T, MPI_COMM_WORLD, &int64MpiPackSize);
        int unsignedIntMpiPackSize;
        MPI_Pack_size(1, MPI_UNSIGNED, MPI_COMM_WORLD, &unsignedIntMpiPackSize);

        // store number of messages;
        int messageSize = unsignedIntMpiPackSize;
        // store 1 int64 per message for flag
        messageSize += num_messages*int64MpiPackSize;
        // store 2 unsigned ints per message for length of tag and length of text
        messageSize += num_messages*2*unsignedIntMpiPackSize;

        for (const auto lm : local_deferredlogger.messages_) {
            int stringMpiPackSize;
            MPI_Pack_size(lm.tag.size(), MPI_CHAR, MPI_COMM_WORLD, &stringMpiPackSize);
            messageSize += stringMpiPackSize;
            MPI_Pack_size(lm.text.size(), MPI_CHAR, MPI_COMM_WORLD, &stringMpiPackSize);
            messageSize += stringMpiPackSize;
        }

        // Pack local messages.
        std::vector<char> buffer(messageSize);

        int offset = 0;
        packMessages(local_deferredlogger.messages_, buffer, offset);
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
        Ewoms::DeferredLogger global_deferredlogger;
        global_deferredlogger.messages_ = unpackMessages(recvBuffer, displ);
        return global_deferredlogger;
    }

} // namespace Ewoms

#else // HAVE_MPI

namespace Ewoms
{
    Ewoms::DeferredLogger gatherDeferredLogger(const Ewoms::DeferredLogger& local_deferredlogger)
    {
        return local_deferredlogger;
    }
} // namespace Ewoms

#endif // HAVE_MPI
