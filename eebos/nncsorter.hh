// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
#ifndef EWOMS_EEBOS_NNCSORTER_HH
#define EWOMS_EEBOS_NNCSORTER_HH

#include <ewoms/eclio/parser/eclipsestate/grid/nnc.hh>
#include <ewoms/eclio/opmlog/opmlog.hh>

#include <ewoms/eclio/opmlog/opmlog.hh>

#include <vector>
#include <algorithm>
#include <sstream>

namespace Ewoms
{
/// \brief Scale NNC data wit informtion form EDITNNC and sort it.
/// \param nncData The NNC data as provided by the deck.
/// \param editnncData The EDITNNC data as provided by the deck.
/// \return A lexicographically sorted vector of the scaled NNC data.
///         For each entry entry.cell1<entry.cell2 will hold for convenience.
inline std::vector<Ewoms::NNCdata> sortNncAndApplyEditnnc(const std::vector<Ewoms::NNCdata>& nncDataIn,
                                                          std::vector<Ewoms::NNCdata> editnncData,
                                                          int numPrintMax = 20)
{
    int numPrinted = 0;

    auto nncLess =
        [](const Ewoms::NNCdata& d1, const Ewoms::NNCdata& d2) {
            return
                (d1.cell1 < d2.cell1)
                || (d1.cell1 == d2.cell1 && d1.cell2 < d2.cell2);
        };

    auto makeCell1LessCell2 =
        [](const Ewoms::NNCdata& entry) {
            if ( entry.cell2 < entry.cell1)
                return Ewoms::NNCdata(entry.cell2, entry.cell1, entry.trans);
            else
                return entry;
        };

    // We need to make sure that for each entry cell1<=cell2 holds. Otherwise sorting
    // will not make the search more accurate if the engineer chooses to define NNCs
    // differently.
    std::vector<Ewoms::NNCdata> nncData(nncDataIn);
    std::transform(nncData.begin(), nncData.end(), nncData.begin(), makeCell1LessCell2);
    std::transform(editnncData.begin(), editnncData.end(), editnncData.begin(), makeCell1LessCell2);
    std::sort(nncData.begin(), nncData.end(), nncLess);
    auto candidate = nncData.begin();

    for (const auto& edit: editnncData) {
        auto printNncWarning =
            [&numPrinted, numPrintMax](int c1, int c2) {
                if (numPrinted < numPrintMax) {
                    std::ostringstream sstr;
                    sstr << "Cannot edit nnc from " << c1 << " to " << c2
                         << " as it does not exist";
                    Ewoms::OpmLog::warning(sstr.str());
                    ++numPrinted;
                }
            };
        if (candidate == nncData.end()) {
            // no more NNCs left
            printNncWarning(edit.cell1, edit.cell2);
            continue;
        }
        if (candidate->cell1 != edit.cell1 || candidate->cell2 != edit.cell2) {
            candidate = std::lower_bound(nncData.begin(), nncData.end(), Ewoms::NNCdata(edit.cell1, edit.cell2, 0), nncLess);
            if (candidate == nncData.end()) {
                // no more NNCs left
                printNncWarning(edit.cell1, edit.cell2);
                continue;
            }
        }
        auto firstCandidate = candidate;
        while (candidate != nncData.end()
               && candidate->cell1 == edit.cell1
               && candidate->cell2 == edit.cell2)
        {
            candidate->trans *= edit.trans;
            ++candidate;
        }
        // start with first match in next iteration to catch case where next
        // EDITNNC is for same pair.
        candidate = firstCandidate;
    }
    return nncData;
}
}
#endif
