# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

cmake_minimum_required(VERSION 3.10)

# set up project and specify the minimum cmake version
project("ewoms-eclsimulators" C CXX)

option(BUILD_EFLOW "Build the production oriented eflow simulator?" OFF)
option(BUILD_EEBOS "Build the research-oriented eebos simulator?" ON)
option(BUILD_EEBOS_EXTENSIONS "Build the variants for various extensions of eebos by default?" ON)
option(BUILD_EEBOS_DEBUG_EXTENSIONS "Build the debugging variants of eebos by default?" OFF)

# find the build system (i.e., ewoms-common) and set cmake's module path
find_package(ewoms-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${ewoms-common_MODULE_PATH})

# include the eWoms cmake macros
include(EwomsMacros)

# do most of the book-keeping required
ewoms_project()

# find the packages needed to compile the unit tests
find_package(Boost OPTIONAL_COMPONENTS unit_test_framework REQUIRED)

# we want all features detected by the build system to be enabled,
# thank you!
dune_enable_all_packages()

ewoms_recursive_add_library("ewomseclsimulators" "ewoms")

# recursively mark all header files beneath the "ewoms" directory for
# installation.
ewoms_recusive_export_all_headers("ewoms")

######
# eebos and its variants
######
# blackoil eebos
ewoms_add_application(eebos
                    DEFAULT_ENABLE_IF ${BUILD_EEBOS}
                    SOURCES eebos/eebos_main.cc $<TARGET_OBJECTS:eebos_libblackoil>)

# eebos variants
set(MEEBOS_TARGETS "")
foreach(EXT blackoil solvent polymer foam gasoil oilwater oilwaterpolymer thermal brine)
  add_library(eebos_lib${EXT} OBJECT EXCLUDE_FROM_ALL eebos/eebos_${EXT}.cc)
  target_include_directories(eebos_lib${EXT} PRIVATE ${ecl_INCLUDE_DIRS})
  list(APPEND MEEBOS_TARGETS $<TARGET_OBJECTS:eebos_lib${EXT}>)
endforeach()

# add specialized eebos simulators with some black-oil extension
foreach(EXT solvent polymer foam gasoil oilwater oilwaterpolymer thermal brine)
  ewoms_add_application(eebos_${EXT}
    DEFAULT_ENABLE_IF ${BUILD_EEBOS_EXTENSIONS}
    SOURCES eebos/eebos_${EXT}_main.cc $<TARGET_OBJECTS:eebos_lib${EXT}>)
endforeach()

######
# meebos: multiplexed eebos with most of the above model extensions
######
ewoms_add_application(meebos
                    DEFAULT_ENABLE_IF ${BUILD_EEBOS_EXTENSIONS}
                    SOURCES eebos/meebos_main.cc ${MEEBOS_TARGETS})

######
# eebos variants only for debugging
######

# variant that uses the default ewoms linear solver and the toy ECL well model
ewoms_add_application(eebos_plain
                    DEFAULT_ENABLE_IF ${BUILD_EEBOS_DEBUG_EXTENSIONS}
                    SOURCES eebos/eebos_plain.cc
                    EXE_NAME eebos_plain)

# variant that uses an AMG-based linear solver
ewoms_add_application(eebos_amg
                    DEFAULT_ENABLE_IF ${BUILD_EEBOS_DEBUG_EXTENSIONS}
                    SOURCES eebos/eebos_blackoil_amg.cc
                    EXE_NAME eebos_amg)

# variant that uses different phase and component indices than the others
ewoms_add_application(eebos_altidx
                    DEFAULT_ENABLE_IF ${BUILD_EEBOS_DEBUG_EXTENSIONS}
                    SOURCES eebos/eebos_altidx.cc
                    EXE_NAME eebos_altidx)

# target to seletively build meebos, eebos, and all eebos variants
add_custom_target(all_eebos)
add_dependencies(all_eebos "eebos")
add_dependencies(all_eebos "meebos")
add_dependencies(all_eebos "eebos_solvent")
add_dependencies(all_eebos "eebos_polymer")
add_dependencies(all_eebos "eebos_foam")
add_dependencies(all_eebos "eebos_thermal")
add_dependencies(all_eebos "eebos_oilwater")
add_dependencies(all_eebos "eebos_gasoil")
add_dependencies(all_eebos "eebos_brine")

######
# eflow (mainly for easier comparison with OPM's 'flow')
######

ewoms_add_application(eflow
  DEFAULT_ENABLE_IF ${BUILD_EFLOW}
  SOURCES
   eflow/eflow.cc
   eflow/eflow_blackoil.cc
   eflow/eflow_gasoil.cc
   eflow/eflow_oilwater.cc
   eflow/eflow_foam.cc
   eflow/eflow_polymer.cc
   eflow/eflow_solvent.cc
   eflow/eflow_energy.cc
   eflow/eflow_brine.cc
   eflow/eflow_oilwater_polymer.cc
   eflow/eflow_oilwater_polymer_injectivity.cc
   )

ewoms_add_application(eflow_blackoil_dunecpr
  DEFAULT_ENABLE_IF ${BUILD_EFLOW}
  SOURCES
  eflow/eflow_blackoil_dunecpr.cc
  )

if(TARGET eflow)
  add_dependencies(all_tests "eflow")
  add_dependencies(test-suite "eflow")
  add_dependencies(build_tests "eflow")
endif()

ewoms_recusive_copy_testdata_to_builddir("tests/*.data" "tests/*.DATA" "tests/VFP*")

ewoms_add_test(test_equil CONDITION Boost_UNIT_TEST_FRAMEWORK_FOUND LIBRARIES "${Boost_LIBRARIES}")
ewoms_add_test(test_ecl_output CONDITION Boost_UNIT_TEST_FRAMEWORK_FOUND LIBRARIES "${Boost_LIBRARIES}")
ewoms_add_test(test_blackoil_amg CONDITION Boost_UNIT_TEST_FRAMEWORK_FOUND LIBRARIES "${Boost_LIBRARIES}")
ewoms_add_test(test_convergencereport CONDITION Boost_UNIT_TEST_FRAMEWORK_FOUND LIBRARIES "${Boost_LIBRARIES}")
ewoms_add_test(test_flexiblesolver CONDITION Boost_UNIT_TEST_FRAMEWORK_FOUND LIBRARIES "${Boost_LIBRARIES}")
ewoms_add_test(test_preconditionerfactory CONDITION Boost_UNIT_TEST_FRAMEWORK_FOUND LIBRARIES "${Boost_LIBRARIES}")
ewoms_add_test(test_graphcoloring CONDITION Boost_UNIT_TEST_FRAMEWORK_FOUND LIBRARIES "${Boost_LIBRARIES}")
ewoms_add_test(test_vfpproperties CONDITION Boost_UNIT_TEST_FRAMEWORK_FOUND LIBRARIES "${Boost_LIBRARIES}")
ewoms_add_test(test_milu CONDITION Boost_UNIT_TEST_FRAMEWORK_FOUND LIBRARIES "${Boost_LIBRARIES}")
ewoms_add_test(test_multmatrixtransposed CONDITION Boost_UNIT_TEST_FRAMEWORK_FOUND LIBRARIES "${Boost_LIBRARIES}")
ewoms_add_test(test_nncsorter CONDITION Boost_UNIT_TEST_FRAMEWORK_FOUND LIBRARIES "${Boost_LIBRARIES}")
ewoms_add_test(test_wellmodel CONDITION Boost_UNIT_TEST_FRAMEWORK_FOUND LIBRARIES "${Boost_LIBRARIES}")
ewoms_add_test(test_deferredlogger CONDITION Boost_UNIT_TEST_FRAMEWORK_FOUND LIBRARIES "${Boost_LIBRARIES}")
ewoms_add_test(test_timer CONDITION Boost_UNIT_TEST_FRAMEWORK_FOUND LIBRARIES "${Boost_LIBRARIES}")
ewoms_add_test(test_invert CONDITION Boost_UNIT_TEST_FRAMEWORK_FOUND LIBRARIES "${Boost_LIBRARIES}")
ewoms_add_test(test_stoppedwells CONDITION Boost_UNIT_TEST_FRAMEWORK_FOUND LIBRARIES "${Boost_LIBRARIES}")
ewoms_add_test(test_relpermdiagnostics CONDITION Boost_UNIT_TEST_FRAMEWORK_FOUND LIBRARIES "${Boost_LIBRARIES}")
ewoms_add_test(test_norne_pvt CONDITION Boost_UNIT_TEST_FRAMEWORK_FOUND LIBRARIES "${Boost_LIBRARIES}")
ewoms_add_test(test_wellstatefullyimplicitblackoil CONDITION Boost_UNIT_TEST_FRAMEWORK_FOUND LIBRARIES "${Boost_LIBRARIES}")
ewoms_add_test(test_parallelistlinformation CONDITION MPI_FOUND AND Boost_UNIT_TEST_FRAMEWORK_FOUND LIBRARIES "${Boost_LIBRARIES}")

# deal with the data required by the unit tests
ewoms_recusive_copy_testdata_to_builddir("tests/*.txt" "tests/VFPPROD2" "tests/*.json" "tests/*.data" "tests/*.DATA")

# finalize the project, e.g. generate the config.h etc.
finalize_ewoms_project()
