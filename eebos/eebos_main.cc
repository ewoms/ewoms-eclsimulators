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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \brief The main file of eebos, the general-purpose black-oil simulator for ECL decks.
 *
 * eebos first parses the input data set before selecting the variant which should run it
 * so it does not require the user to select the correct binary. On the flipside, only
 * the most common options can be supported by the simulator to avoid a combinatorial
 * explosion of combinations and even then compiling the simulator takes mighty
 * long. (I.e., consider using the a specialized binary for development.) The options
 * available by default are standard blackoil, twophase oil-water, twophase gas-oil,
 * solvent, polymer, foam and thermal.
 */
#include "config.h"

#include "eebos_blackoil.hh"
#include "eebos_oilwater.hh"
#include "eebos_gasoil.hh"
// TODO (?): #include "eebos_watergas.hh"
#include "eebos_thermal.hh"
#include "eebos_solvent.hh"
#include "eebos_polymer.hh"
#include "eebos_foam.hh"

#include <ewoms/common/propertysystem.hh>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>

#include <ewoms/eclio/parser/parser.hh>
#include <ewoms/eclio/parser/parsecontext.hh>
#include <ewoms/eclio/parser/errorguard.hh>
#include <ewoms/eclio/parser/deck/deck.hh>

#include <string>
#include <memory>

int main(int argc, char **argv)
{
    Dune::Timer externalSetupTimer;
    externalSetupTimer.start();

    if (!Ewoms::eebosBlackOilDeckFileNameIsSet(argc, argv))
        // no deck was specified, e.g., --help. use the black oil variant to figure out
        // what exactly should be done
        return Ewoms::eebosBlackOilMain(argc, argv);

    std::string deckFileName =
        Ewoms::eebosBlackOilGetDeckFileName(argc, argv);

    std::unique_ptr<Ewoms::ParseContext> parseContext
        = Ewoms::eebosBlackOilCreateParseContext(argc, argv);
    std::unique_ptr<Ewoms::ErrorGuard> errorGuard(new Ewoms::ErrorGuard);

    // deal with parallel runs
    int myRank = Dune::MPIHelper::instance(argc, argv).rank();

    Ewoms::Parser parser;
    // parse the deck file
    if (myRank == 0)
        std::cout << "Parsing deck file \"" << deckFileName << "\"" << std::endl;
    std::unique_ptr<Ewoms::Deck> deck(new Ewoms::Deck(parser.parseFile(deckFileName, *parseContext, *errorGuard)));

    // TODO: check which variant ought to be used
    bool waterActive = deck->hasKeyword("WATER");
    bool gasActive = deck->hasKeyword("GAS");
    bool oilActive = deck->hasKeyword("OIL");
    bool solventActive = deck->hasKeyword("SOLVENT");
    bool polymerActive = deck->hasKeyword("POLYMER");
    bool foamActive = deck->hasKeyword("FOAM");
    bool thermalActive = deck->hasKeyword("THERMAL") || deck->hasKeyword("TEMP");

    std::stringstream notSupportedErrorStream;
    notSupportedErrorStream << "deck not supported by eebos, you might want to use a specialized binary. Active options:\n"
                            << "   water: " << waterActive << "\n"
                            << "   gas: " << gasActive << "\n"
                            << "   oil: " << oilActive << "\n"
                            << "   solvent: " << solventActive << "\n"
                            << "   polymer: " << polymerActive << "\n"
                            << "   foam: " << foamActive << "\n"
                            << "   thermal/temperature: " << thermalActive << "\n";

    int numBlackOilPhases = (waterActive?1:0) + (gasActive?1:0) + (oilActive?1:0);
    if (numBlackOilPhases == 0) {
        notSupportedErrorStream << "\n"
                                << "no black-oil phase (water, gas or oil) specified.\n";
        std::cerr << notSupportedErrorStream.str() << std::endl;
        std::abort();
    }
    else if (numBlackOilPhases == 1) {
        notSupportedErrorStream << "\n"
                                << "single-phase simulations are unsupported\n";
        std::cerr << notSupportedErrorStream.str() << std::endl;
        std::abort();
    }
    else if (numBlackOilPhases == 2) {
        if (solventActive) {
            notSupportedErrorStream << "\n"
                                    << "combining twophase and solvent is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (polymerActive) {
            notSupportedErrorStream << "\n"
                                    << "combining twophase and polymer is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (foamActive) {
            notSupportedErrorStream << "\n"
                                    << "combining twophase and foam is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (thermalActive) {
            notSupportedErrorStream << "\n"
                                    << "combining twophase and energy conservation is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (oilActive && waterActive) {
            if (myRank == 0)
                std::cout << "Using oil-water mode" << std::endl;
            Ewoms::eebosOilWaterSetDeck(deck.get(),
                                     parseContext.get(),
                                     errorGuard.get(),
                                     externalSetupTimer.elapsed());
            return Ewoms::eebosOilWaterMain(argc, argv);
        }
        else if (oilActive && gasActive) {
            // run eebos_gasoil
            if (myRank == 0)
                std::cout << "Using gas-oil mode" << std::endl;
            Ewoms::eebosGasOilSetDeck(deck.get(),
                                   parseContext.get(),
                                   errorGuard.get(),
                                   externalSetupTimer.elapsed());
            return Ewoms::eebosGasOilMain(argc, argv);
        }
        else if (waterActive && gasActive) {
            notSupportedErrorStream << "\n"
                                    << "water-gas simulations are currently unsupported\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }
    }
    else if (foamActive) {
        if (solventActive) {
            notSupportedErrorStream << "\n"
                                    << "combining foam and solvent is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (polymerActive) {
            notSupportedErrorStream << "\n"
                                    << "combining foam and polymer is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (thermalActive) {
            notSupportedErrorStream << "\n"
                                    << "combining foam and and energy conservation is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        // run eebos_foam
        if (myRank == 0)
            std::cout << "Using foam mode" << std::endl;
        Ewoms::eebosFoamSetDeck(deck.get(),
                             parseContext.get(),
                             errorGuard.get(),
                             externalSetupTimer.elapsed());
        return Ewoms::eebosFoamMain(argc, argv);
    }
    else if (polymerActive) {
        if (solventActive) {
            notSupportedErrorStream << "\n"
                                    << "combining polymer and solvent is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (foamActive) {
            notSupportedErrorStream << "\n"
                                    << "combining polymer and foam is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (thermalActive) {
            notSupportedErrorStream << "\n"
                                    << "combining polymer and and energy conservation is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        // run eebos_polymer
        if (myRank == 0)
            std::cout << "Using polymer mode" << std::endl;
        Ewoms::eebosPolymerSetDeck(deck.get(),
                                parseContext.get(),
                                errorGuard.get(),
                                externalSetupTimer.elapsed());
        return Ewoms::eebosPolymerMain(argc, argv);
    }
    else if (solventActive) {
        if (polymerActive) {
            notSupportedErrorStream << "\n"
                                    << "combining solvent and polymer is not supported\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (foamActive) {
            notSupportedErrorStream << "\n"
                                    << "combining solvent and foam is not supported\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (thermalActive) {
            notSupportedErrorStream << "\n"
                                    << "combining solvent and and energy conservation is not supported\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        // run eebos_solvent
        if (myRank == 0)
            std::cout << "Using solvent mode" << std::endl;
        Ewoms::eebosSolventSetDeck(deck.get(),
                                parseContext.get(),
                                errorGuard.get(),
                                externalSetupTimer.elapsed());
        return Ewoms::eebosSolventMain(argc, argv);
    }
    else if (thermalActive) {
        if (solventActive) {
            notSupportedErrorStream << "\n"
                                    << "combining thermal and solvent is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (polymerActive) {
            notSupportedErrorStream << "\n"
                                    << "combining thermal and polymer is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        if (foamActive) {
            notSupportedErrorStream << "\n"
                                    << "combining thermal and foam is not supported by the multiplexed simulator\n";
            std::cerr << notSupportedErrorStream.str() << std::endl;
            std::abort();
        }

        // run eebos_thermal
        if (myRank == 0)
            std::cout << "Using thermal mode" << std::endl;
        Ewoms::eebosThermalSetDeck(deck.get(),
                                parseContext.get(),
                                errorGuard.get(),
                                externalSetupTimer.elapsed());
        return Ewoms::eebosThermalMain(argc, argv);
    }
    else {
        if (myRank == 0)
            std::cout << "Using blackoil mode" << std::endl;
        Ewoms::eebosBlackOilSetDeck(deck.get(),
                                 parseContext.get(),
                                 errorGuard.get(),
                                 externalSetupTimer.elapsed());
        return Ewoms::eebosBlackOilMain(argc, argv);
    }

    if (myRank == 0)
        // this is supposed to be unreachable. this should not happen!
        std::cerr << "Oops: something went wrong when deciding which simulator ought to be used" << std::endl;
    std::abort();

    return 0;
}
