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
#ifndef EWOMS_MISSINGFEATURES_HH
#define EWOMS_MISSINGFEATURES_HH

#include <ewoms/eclio/parser/deck/deck.hh>

#include <string>
#include <map>

namespace Ewoms {

class ErrorGuard;
class ParseContext;

namespace MissingFeatures {

    template <typename T>
    struct PartiallySupported {
        std::string item;
        T item_value;
    };

    template <typename Keyword, typename Item, typename T>
    void addSupported(std::multimap<std::string, PartiallySupported<T> >& map, T itemValue);

    template <typename T>
    void checkOptions(const DeckKeyword& keyword, std::multimap<std::string , PartiallySupported<T> >& map, const ParseContext& parseContext, ErrorGuard& errorGuard);

    void checkKeywords(const Deck& deck, const ParseContext& parseContext, ErrorGuard& errorGuard);

    template<typename T>
    void checkKeywords(const Deck& deck, const ParseContext& parseContext, T&& errorGuard);

    void checkKeywords(const Deck& deck);
}

}

#endif // EWOMS_MISSINGFEATURES_HH
