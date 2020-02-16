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

#ifndef EWOMS_PARALLELFILEMERGER_HH
#define EWOMS_PARALLELFILEMERGER_HH

#include <memory>
#include <iostream>

#include <ewoms/common/filesystem.hh>
#include <regex>

namespace Ewoms
{
namespace detail
{

namespace fs = Ewoms::filesystem;

/// \brief A functor that merges multiple files of a parallel run to one file.
///
/// Without care multiple processes might log messages in a parallel run.
/// Non-root processes will do that to seperate files
/// <basename>.<rank>.<extension. This functor will append those file
/// to usual ones and delete the other files.
class ParallelFileMerger
{
public:
    /// \brief Constructor
    /// \param output_dir The output directory to use for reading/Writing.
    /// \param deckanme The name of the deck.
    ParallelFileMerger(const fs::path& output_dir,
                       const std::string& deckname,
                       bool show_fallout = false)
        : debugFileRegex_(deckname+"\\.\\d+\\.DBG"),
          logFileRegex_(deckname+"\\.\\d+\\.PRT"),
          fileWarningRegex_(deckname+"\\.(\\d+)\\.[^.]+"),
          show_fallout_(show_fallout)
    {
        if ( show_fallout_ )
        {
            auto debugPath = output_dir;
            debugPath /= (deckname + ".DBG");
            debugStream_.reset(new std::ofstream(debugPath.string(),
                                                 std::ofstream::app));
            auto logPath = output_dir;
            logPath /= ( deckname + ".PRT");
            logStream_.reset(new std::ofstream(logPath.string(),
                                               std::ofstream::app));
        }
    }

    void operator()(const fs::path& file)
    {
        std::smatch matches;
        std::string filename = file.filename().native();

        if ( std::regex_match(filename, matches, fileWarningRegex_) )
        {
            std::string rank = std::regex_replace(filename, fileWarningRegex_, "\\1");

            if( std::regex_match(filename, logFileRegex_) )
            {
                if ( show_fallout_ ){
                    appendFile(*logStream_, file, rank);
                }else{
                    fs::remove(file);
                }
            }
            else
            {
                if (std::regex_match(filename, debugFileRegex_)  )
                {
                    if ( show_fallout_ ){
                        appendFile(*debugStream_, file, rank);
                    }else{
                        fs::remove(file);
                    }
                }
                else
                {
                    if ( show_fallout_ ){
                        std::cerr << "WARNING: Unrecognized file with name "
                                  << filename
                                  << " that might stem from a  parallel run."
                                  << std::endl;
                    }
                }
            }
        }
    }
private:
    /// \brief Append contents of a file to a stream
    /// \brief of The output stream to use.
    /// \brief file The file whose content to append.
    /// \brief rank The rank that wrote the file.
    void appendFile(std::ofstream& of, const fs::path& file, const std::string& rank)
    {
        if( fs::file_size(file) )
        {
            std::cerr << "WARNING: There has been logging to file "
                      << file.string() <<" by process "
                      << rank << std::endl;

            std::ifstream in(file);
            of<<std::endl<< std::endl;
            of<<"=======================================================";
            of<<std::endl<<std::endl;
            of << " Output written by rank " << rank << " to file " << file.string();
            of << ":" << std::endl << std::endl;
            of << in.rdbuf() << std::endl << std::endl;
            of << "======================== end output =====================";
            of << std::endl;
            in.close();
        }
        fs::remove(file);
    }

    /// \brief Regex to capture *.DBG
    std::regex debugFileRegex_;
    /// \brief Regex to capture  *.PRT
    std::regex logFileRegex_;
    /// \brief Regex to capture  CASENAME.[0-9]+.[A-Z]+
    std::regex fileWarningRegex_;
    /// \brief Stream to *.DBG file
    std::unique_ptr<std::ofstream> debugStream_;
    /// \brief Stream to *.PRT file
    std::unique_ptr<std::ofstream> logStream_;
    /// \brief Whether to show any logging fallout
    bool show_fallout_;
};
} // end namespace detail
} // end namespace Ewoms
#endif // end header guard
