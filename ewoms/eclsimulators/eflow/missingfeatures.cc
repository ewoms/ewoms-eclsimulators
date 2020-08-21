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

#include <ewoms/eclsimulators/eflow/missingfeatures.hh>

#include <ewoms/eclio/opmlog/opmlog.hh>

#include <ewoms/eclio/parser/deck/deckkeyword.hh>
#include <ewoms/eclio/parser/parser.hh>
#include <ewoms/eclio/parser/parsecontext.hh>
#include <ewoms/eclio/parser/errorguard.hh>
#include <ewoms/eclio/parser/parserkeywords/c.hh>
#include <ewoms/eclio/parser/parserkeywords/e.hh>
#include <ewoms/eclio/parser/parserkeywords/p.hh>

#include <unordered_set>
#include <string>
#include <map>
#include <type_traits>

namespace Ewoms {

namespace MissingFeatures {

    template <typename Keyword, typename Item, typename T>
    void addSupported(std::multimap<std::string, PartiallySupported<T> >& map, T itemValue)
    {
        std::pair<std::string,PartiallySupported<T> > pair({Keyword::keywordName, PartiallySupported<T>{Item::itemName , itemValue}});
        map.insert(pair);
    }

template <class T>
std::string optionToString_(std::enable_if_t<std::is_arithmetic<T>::value, T>& val)
{ return std::to_string(val); }

template <class T>
std::string optionToString_(std::enable_if_t<!std::is_arithmetic<T>::value, T>& val)
{ return val; }

    template <typename T>
    void checkOptions(const DeckKeyword& keyword, std::multimap<std::string , PartiallySupported<T> >& map, const ParseContext& parseContext, ErrorGuard& errorGuard)
    {
        // check for partially supported keywords.
        typename std::multimap<std::string, PartiallySupported<T> >::iterator it, itlow, itup;
        itlow = map.lower_bound(keyword.name());
        itup  = map.upper_bound(keyword.name());
        for (it = itlow; it != itup; ++it) {
            const auto& record = keyword.getRecord(0);
            if (record.getItem(it->second.item).template get<T>(0) != it->second.item_value) {
                const auto& location = keyword.location();
                std::string val = optionToString_<decltype(it->second.item_value)>(it->second.item_value);

                std::string msg = "For keyword '" + it->first + "' only value " + val
                    + " in item " + it->second.item + " is supported by eflow.\n"
                    + "In file " + location.filename + ", line " + std::to_string(location.lineno) + "\n";
                parseContext.handleError(ParseContext::SIMULATOR_KEYWORD_ITEM_NOT_SUPPORTED, msg, errorGuard);
            }
        }
    }

    template <typename T>
    void checkKeywords(const Deck& deck, const ParseContext& parseContext, T&& errorGuard) {
        checkKeywords(deck, parseContext, errorGuard);
    }

    void checkKeywords(const Deck& deck) {
        checkKeywords(deck, ParseContext(), ErrorGuard());
    }

    void checkKeywords(const Deck& deck, const ParseContext& parseContext, ErrorGuard& errorGuard)
    {
        // These keywords are supported by ewoms-eclio, but are not supported
        // by eflow. For some of them, only part of the options are supported.
        // The list is used to output messages only.
        std::unordered_set<std::string> unsupported_keywords = {
            "ACTION",
            "ACTIONG",
            "ACTIONR",
            "ACTIONS",
            "ACTIONW",
            "ACTPARAM",
            "ADSALNOD",
            "ADDZCORN",
            "ADSORP",
            "AITS",
            "AITSOFF",
            "ALKADS",
            "ALKALINE",
            "ALKROCK",
            "API",
            "ALPOLADS",
            "ALSURFAD",
            "ALSURFST",
            "AMALGAM",
            "APIGROUP",
            "APILIM",
            "APIVD",
            "AQANCONL",
            "AQANNC",
            "AQANTRC",
            "AQUALIST",
            "AQUCHGAS",
            "AQUCHWAT",
            "AQUCON",
            "AQUCWFAC",
            "AQUFET",
            "AQUFLUX",
            "AQUNNC",
            "AQUNUM",
            "AUTOCOAR",
            "AUTOREF",
            "BIGMODEL",
            "BDENSITY",
            "BGGI",
            "BOGI",
            "BOUNDARY",
            "BPARA",
            "BPIDIMS",
            "BRANPROP",
            "BTOBALFA",
            "BTOBALFV",
            "CALTRAC",
            "CARFIN",
            "CART",
            "CBMOPTS",
            "CECON",
            "CECONT",
            "COAL",
            "COALADS",
            "COALNUM",
            "COALPP",
            "COARSEN",
            "COLLAPSE",
            "COLUMNS",
            "CBMOPTS",
            "COMPDATX",
            "COMPDATL",
            "COMPDATM",
            "COMPDATL",
            "COMPIMB",
            "COMPFLSH",
            "COMPINJK",
            "COMPLMPL",
            "COMPOFF",
            "COMPRIV",
            "COMPRP",
            "COMPRPL",
            "COMPSEGL",
            "COMPVE",
            "COMPVEL",
            "CPIFACT",
            "CPIFACTL",
            "CPR",
            "CSKIN",
            "CONNECTION",
            "CONNECTION_PROBE",
            "COORDSYS",
            "COPYBOX",
            "COPYREG",
            "CRITPERM",
            "DATUMR",
            "DATUMRX",
            "DCQDEFN",
            "DEBUG",
            "DELAYACT",
            "DEPTHTAB",
            "DIAGDISP",
            "DIFF",
            "DIFFC",
            "DIFFCOAL",
            "DIFFDP",
            "DIFFMMF",
            "DIFFMR",
            "DIFFMR-",
            "DIFFMTHT",
            "DIFFMTH-",
            "DIFFMX",
            "DIFFMX-",
            "DIFFMY",
            "DIFFMY-",
            "DIFFMZ",
            "DIFFMZ-",
            "DIFFR",
            "DIFFTHT",
            "DIFFUSE"
            "DIFFX",
            "DIFFY",
            "DIFFZ",
            "DIMPES",
            "DIMPLICT",
            "DISPDIMS",
            "DISPERSE",
            "DOMAINS",
            "DPGRID",
            "DPKRMOD",
            "DPNUM",
            "DR",
            "DRV",
            "DRILPRI",
            "DSPDEINT",
            "DUALPERM",
            "DUALPORO",
            "DUMPCUPL",
            "DUMPFLUX",
            "DYNAMICR",
            "DYNRDIMS",
            "DYNAMICR",
            "DZMATRIX",
            "DZMTRX",
            "DZMTRXV",
            "DZNET",
            "ECHO",
            "ECLMC",
            "EDITNNCR",
            "EHYSTRR",
            "ENDDYN",
            "ENDFIN",
            "ENDNUM",
            "ENDSKIP",
            "ENKRVD",
            "ENKSRVD",
            "ENPCVD",
            "ENPTVD",
            "ENSPCVD",
            "EPSDBGS",
            "EPSDEBUG",
            "EQLZCORN",
            "EQUALREG",
            "ESSNODE",
            "EXCAVATE",
            "EXCEL",
            "EXTFIN",
            "EXTHOST",
            "EXTRAPMS",
            "EXTREPGL",
            "FBHPDEF",
            "FHERCHBL",
            "FRICTION",
            "FIPSEP",
            "FLUXREG",
            "FLUXTYPE",
            "FMTHMD",
            "FOAMDCYO",
            "FOAMDCYW",
            "FOAMFCN",
            "FOAMFRM",
            "FOAMFSO"
            "FOAMFST",
            "FOAMFSW",
            "FOAMMOBP",
            "FOAMMOBS",
            "FORMFEED",
            "FULLIMP",
            "GEGONT",
            "GETDATA",
            "GASBEGIN",
            "GASCONC",
            "GASEND",
            "GASFCOMP",
            "GASFDECR",
            "GASFDELC",
            "GASFIELD",
            "GASFTARG",
            "GASMONTH",
            "GASPERIO",
            "GASSATC",
            "GASYEAR",
            "GCALECON",
            "GCONCAL",
            "GCONENG",
            "GCONPRI",
            "GCONSALE",
            "GCONSUMP",
            "GCONTOL",
            "GCUTBACK",
            "GCUTBACT",
            "GCVD",
            "GDCQ",
            "GDCQECON",
            "GDIMS",
            "GDORIENT",
            "GDRILPOT",
            "GECON",
            "GECONT",
            "GETGLOB",
            "GI",
            "GIALL",
            "GIMODEL",
            "GINODE",
            "GLIFTLIM",
            "GLIFTOPT",
            "GNETDP",
            "GNETINJE",
            "GNETPUMP",
            "GPMAINT",
            "GRADGRUP",
            "GRADRESV",
            "GRADRFT",
            "GRADWELL",
            "GRAVCONS",
            "GRAVDR",
            "GRAVDRB",
            "GRAVDRM",
            "GRDREACH",
            "GRIDUNIT",
            "GRUPMAST",
            "GRUPNET",
            "GRUPRIG",
            "GRUPSLAV",
            "GRUPTARG",
            "GSATINJE",
            "GSEPCOND",
            "GSSCPTST",
            "GSWINGF",
            "GTADD",
            "GTMULT",
            "GUIDECAL",
            "GSATPROD",
            "GUIDERAT",
            "GUPFREQ",
            "GWRTWCV",
            "HALFTRAN",
            "HAxxxxxx",
            "HBNUM",
            "HDISP",
            "HMAQUCT",
            "HMAQUFET",
            "HMAQUNUM",
            "HMDIMS",
            "HMFAULTS",
            "HMMLAQUN",
            "HMMLCTAQ",
            "HMMLFTAQ",
            "HMMLTWCN",
            "HMMULTxx",
            "HMMULTFT",
            "HMMULTSG",
            "HMPROPS",
            "HMROCK",
            "HMROCKT",
            "HMRREF",
            "HMWELCON",
            "HMWPIMLT",
            "HMxxxxxx",
            "HRFIN",
            "HWKRO",
            "HWKRORG",
            "HWKRORW",
            "HWKRW",
            "HWKRWR",
            "HWPCW",
            "HWSNUM",
            "HWSOGCR",
            "HWSOWCR",
            "HWSWCR",
            "HWSWL",
            "HWSWLPC",
            "HWSWU",
            "HXFIN",
            "HYDRHEAD",
            "HYFIN",
            "HYMOBGDR",
            "HYST",
            "HYSTCHCK",
            "HZFIN",
            "IHOST",
            "IMBNUMMF",
            "IMKRVD",
            "IMPCVD",
            "IMPES",
            "IMPLICIT",
            "IMPORT",
            "IMPTVD",
            "INRAD",
            "IMSPCVD",
            "INSPEC",
            "INTPC",
            "IONROCK",
            "IONXROCK",
            "IONXSURF",
            "ISOLNUM",
            "JFUNCR",
            "KRNUM",
            "KRNUMMF",
            "LANGMPL",
            "LANGMUIR",
            "LANGSOLV",
            "LCUNIT",
            "LGR",
            "LGRCOPY",
            "LGRFREE",
            "LGRLOCK",
            "LGROFF",
            "LGRON",
            "LICENSES",
            "LINCOM",
            "LINKPERM",
            "LIFTOPT",
            "LKRO",
            "LKRORG",
            "LKRORW",
            "LKRW",
            "LKRWR",
            "LOAD",
            "LOWSALT",
            "LPCW",
            "LSALTFNC",
            "LSLTWNUM",
            "LSNUM",
            "LSOGCR",
            "LSOWCR",
            "LSWCR",
            "LSWL",
            "LSWLPC",
            "LSWU",
            "LTOSIGMA",
            "LWKRO",
            "LWKRORG",
            "LWKRORW",
            "LWKRW",
            "LWKRWR",
            "LWPCW",
            "LWSLTNUM",
            "LWSNUM",
            "LWSOGCR",
            "LWSOWCR",
            "LWSWCR",
            "LWSWL",
            "LWSWLPC",
            "LWSWU",
            "LX",
            "LXFIN",
            "LY",
            "LYFIN",
            "LZ",
            "LZFIN",
            "MAPUNITS",
            "MASSFLOW",
            "MATCORR",
            "MEMORY",
            "MESSAGE",
            "MESSOPTS",
            "MESSSRVC",
            "MINNNCT",
            "MINPORV",
            "MLANG",
            "MLANGSLV",
            "MONITOR",
            "MPFANUM",
            "MPFNNC",
            "MSGFILE",
            "MULSGGD",
            "MULSGGDV",
            "MULTOUTS",
            "MULTREAL",
            "MULTREGD",
            "MULTREGH",
            "MULTSIG",
            "MULTSIGV",
            "MULT_XYZ",
            "NARROW",
            "NCONSUMP",
            "NEFAC",
            "NETBALAN",
            "NETCOMPA",
            "NETWORK",
            "NEXT",
            "NEXTSTEP",
            "NEXTSTPL",
            "NINENUM",
            "NINEPOIN",
            "NMATOPTS",
            "NMATRIX",
            "NODPPM",
            "NOECHO",
            "NOHMD",
            "NOHMO",
            "NOHYST",
            "NOWARNEP",
            "NRSOUT",
            "NNEWTF",
            "NOCASC",
            "NODEPROP",
            "NOGGF",
            "NOINSPEC",
            "NLINEARS",
            "NOMONITO",
            "NONNC",
            "NORSSPEC",
            "NOWARN",
            "NSTACK",
            "NUMRES",
            "NWATREM",
            "NXFIN",
            "NYFIN",
            "NZFIN",
            "OFM",
            "OILAPI",
            "OLDTRAN",
            "OLDTRANR",
            "OPERATER",
            "OPERNUM",
            "OPTIONS",
            "OPTIONS",
            "OUTSOL",
            "PARAOPTS",
            "PCG32D",
            "PCW32D",
            "PERMJFUN",
            "PETOPTS",
            "PLYESAL",
            "PLYKRRF",
            "PLYOPTS",
            "PLYRMDEN",
            "PLYROCKM",
            "PLYTRRF",
            "PLYTRRFA",
            "PLYVISCS",
            "PLYVISCT",
            "PLYVSCST",
            "PVZG",
            "PMAX",
            "PRIORITY",
            "PSTEADY",
            "PSWRG",
            "PSWRO",
            "PVCO",
            "PVZG",
            "QDRILL",
            "QDRILL",
            "QHRATING",
            "QMOBIL",
            "PARALLEL",
            "PARTTRAC",
            "PBUB",
            "PCG",
            "PCW",
            "PDEW",
            "PEBI",
            "PECOEFS",
            "PEDIMS",
            "PEGTABX",
            "PEKTABX",
            "PENUM",
            "PERMAVE",
            "PERMFACT",
            "PERMR",
            "PERMTHT",
            "PERMXY",
            "PERMYZ",
            "PERMZX",
            "PETGRID",
            "PICOND",
            "PIMULTAB",
            "PINCHNUM",
            "PINCHOUT",
            "PINCHREG",
            "PINCHXY",
            "PLYADSS",
            "PLYATEMP",
            "PLYCAMAX",
            "PLYDHFLF",
            "PPCWMAX",
            "PRECSALT",
            "PRORDER",
            "PRVD",
            "PVTGW",
            "PVTGWO",
            "RAINFALL",
            "RBEDCONT",
            "RADFIN",
            "RADFIN4",
            "RADIAL",
            "RCMASTS",
            "REACACT",
            "REACHES",
            "READDATA",
            "RESIDNUM",
            "RESVNUMS",
            "RIVDEBUG",
            "RIVRXSEC",
            "RIVERSYS",
            "RIVRDIMS",
            "RIVRPROP",
            "RIVRXSE",
            "RIVSALT",
            "RIVTRACE",
            "ROCKFRAC",
            "ROCKPAMA",
            "ROCKTAB",
            "ROCKTABH",
            "ROCKTABW",
            "ROCKTHSG",
            "ROCKTSIG",
            "ROCKV",
            "RPTCPL",
            "RPTGRIDL",
            "RPTHM",
            "RPTHMG",
            "RPTHMD",
            "RPTHMW",
            "RPTINIT",
            "RPTISOL",
            "RPTPROPS",
            "RPTREGS",
            "RSGI",
            "RSSPE",
            "RSSSPEC",
            "RVCONS",
            "RVCONSTT",
            "RVGI",
            "REFINE",
            "RADFIN4",
            "RHO",
            "RKTRMDIR",
            "ROCKOPTS",
            "ROCKTAB",
            "RPTGRID",
            "RPTONLY",
            "RPTONLYO",
            "RPTPROS",
            "PRTRST",
            "RPTRUNSP",
            "RPTSMRY",
            "RPTSOL",
            "RSCONST",
            "RSCONSTT",
            "RSSPEC",
            "RTEMPA",
            "RWGSALT",
            "SALTPVD",
            "SALTSOL",
            "SAMG",
            "SAVE",
            "SKIP",
            "SKIP100",
            "SKIP300",
            "SUMTHIN",
            "SALT",
            "SALTNODE",
            "SALTREST",
            "SCALELIM",
            "SCDATAB",
            "SCDETAB",
            "SCDPTAB",
            "SCDPTRAC",
            "SCDPDIMS",
            "SCVD",
            "SEPVALS",
            "SFOAM",
            "SGF32D",
            "SIGMA",
            "SIGMAGD",
            "SIGMAGDV",
            "SIGMATH",
            "SIGMAV",
            "SIMULATE",
            "SKRO",
            "SKRORG",
            "SKRORW",
            "SKRW",
            "SKRWR",
            "SLAVES",
            "SMULTX",
            "SMULTY",
            "SMULTZ",
            "SOCRS",
            "SOF32D",
            "SOLVCONC",
            "SOLVDIMS",
            "SOLVDIRS",
            "SOLVFRAC",
            "SOLVNUM",
            "SOLWNUM",
            "SOMGAS",
            "SOMWAT",
            "SSGCR",
            "SSGL",
            "SSOGCR",
            "SSOWCR",
            "SSWCR",
            "SSWL",
            "SSWU",
            "STOG",
            "STOW",
            "STWG",
            "SURF",
            "SURFACT",
            "SURFACTW",
            "SURFADDW",
            "SURFADS",
            "SURFCAPD",
            "SURFESAL",
            "SURFNUM",
            "SURFOPTS",
            "SURFROCK",
            "SURFST",
            "SURFSTES",
            "SURFVISC",
            "SURFWNUM",
            "SWF32D",
            "SWINGFAC"
            "TEMP",
            "TEMPNODE",
            "TEMPTVD",
            "TEMPVD",
            "TIGHTEN",
            "TIGHTENP",
            "TIME",
            "TNUM",
            "TPAMEPS",
            "TPAMEPSS",
            "TRACERKM",
            "TRACERKP",
            "TRACITVD",
            "TRACTVD",
            "TRACER",
            "TRACERS",
            "TRACITVD",
            "TRADS",
            "TRANGL",
            "TRANR",
            "TRANTHT",
            "TRDCY",
            "TRDIF",
            "TRDIS",
            "TRKPF",
            "TRNHD",
            "TRPLPORO",
            "TRROCK",
            "TUNINGDP",
            "TUNINGH",
            "TUNINGL",
            "TUNINGS",
            "TVDP",
            "TZONE",
            "UDT",
            "UDTDIMS",
            "UNCODHMD",
            "UNIFOUTS",
            "UNIFSAVE",
            "USECUPL",
            "USEFLUX",
            "USENOFLO",
            "VAPWAT",
            "VDFLOW",
            "VDFLOWR",
            "VE",
            "VEDEBUG",
            "VEFIN",
            "VEFRAC",
            "VEFRACP",
            "VEFRACPV",
            "VEFRACV",
            "VFPCHK",
            "VFPTABL",
            "VISAGE",
            "VISCD",
            "VISDATES",
            "VISOPTS",
            "WAGHYSTR",
            "WAITBAL",
            "WALKALIN",
            "WALQCALC",
            "WAPI",
            "WARN",
            "WBHGLR",
            "WBOREVOL",
            "WCALCVAL",
            "WBP",
            "WBP4",
            "WBP5",
            "WBP9",
            "WCONINJ",
            "WCONINJP",
            "WCUTBACK",
            "WCUTBACT",
            "WCYCLE",
            "WDFACCOR",
            "WDFAC",
            "WDRILTIM",
            "WDRILPRI",
            "WDRILRES",
            "WECONINJ",
            "WECONT",
            "WEIR",
            "WEIT",
            "WELCNTL",
            "WELDEBUG",
            "WELDRAW",
            "WELEVNT",
            "WELMOVEL"
            "WELOPENL"
            "WELPI",
            "WELPRI",
            "WELSOMIN",
            "WELSPECL",
            "WFRICSEG",
            "WFRICSGL",
            "WFRICTN",
            "WFRICTNL",
            "WGASPROD",
            "WGORPEN",
            "WH2NUM",
            "WH3NUM",
            "WHEDREFD",
            "WHTEMP",
            "WINJMULT",
            "WLIMTOL",
            "WLIFT",
            "WLIFTOPT",
            "WLISTARG",
            "WLISTNAM",
            "WLISTOPT",
            "WNETCTRL",
            "WNETDP",
            "WORKLIM",
            "WORKTHP",
            "WPAVE",
            "WPIMULTL",
            "WPITAB",
            "WPLUG",
            "WPOLYRED",
            "WPOTCALC",
            "WREGROUP",
            "WSCCLEAN",
            "WSCCLENL",
            "WSCTAB",
            "WSEGAICD",
            "WSEGDFIN",
            "WSEGDFMD",
            "WSEGDFPA",
            "WSEGEXSS",
            "WSEGFLIM",
            "WSEGFMOD",
            "WSEGINIT",
            "WSEGITER",
            "WSEGLABY",
            "WSEGLINK",
            "WSEGMULT",
            "WSEGPROP",
            "WSEGPULL",
            "WSEGSEP",
            "WSEGSOLV",
            "WSEGTABL",
            "WSURFACT",
            "WTADD",
            "WTEMPQ",
            "WTHPMAX",
            "WTMULT",
            "WVFPDP",
            "WVFPEXP",
            "WWPAVE",
            "ZIPPY2",
            "ZIPP2OFF"};

        std::multimap<std::string, PartiallySupported<std::string> > string_options;
        std::multimap<std::string, PartiallySupported<int> > int_options;
        addSupported<ParserKeywords::COMPORD, ParserKeywords::COMPORD::ORDER_TYPE, std::string>(string_options , "INPUT");
        addSupported<ParserKeywords::ENDSCALE, ParserKeywords::ENDSCALE::DIRECT, std::string>(string_options, "NODIR");
        addSupported<ParserKeywords::ENDSCALE, ParserKeywords::ENDSCALE::IRREVERS, std::string>(string_options, "REVER");
        addSupported<ParserKeywords::PINCH, ParserKeywords::PINCH::CONTROL_OPTION, std::string>(string_options, "GAP");
        addSupported<ParserKeywords::PINCH, ParserKeywords::PINCH::PINCHOUT_OPTION, std::string>(string_options, "TOPBOT");
        addSupported<ParserKeywords::EHYSTR, ParserKeywords::EHYSTR::relative_perm_hyst, int>(int_options , 0);

        // check deck and keyword for eflow and parser.
        for (size_t idx = 0; idx < deck.size(); ++idx) {
            const auto& keyword = deck.getKeyword(idx);
            const auto& location = keyword.location();
            std::unordered_set<std::string>::const_iterator it;
            it = unsupported_keywords.find(keyword.name());
            if (it != unsupported_keywords.end()) {
                std::string msg = "Keyword '" + keyword.name() + "' is not supported by eflow.\n"
                    + "In file " + location.filename + ", line " + std::to_string(location.lineno) + "\n";
                parseContext.handleError(ParseContext::SIMULATOR_KEYWORD_NOT_SUPPORTED, msg, errorGuard);
            }
            checkOptions<std::string>(keyword, string_options, parseContext, errorGuard);
            checkOptions<int>(keyword, int_options, parseContext, errorGuard);
        }
    }
} // namespace MissingFeatures
} // namespace Ewoms
