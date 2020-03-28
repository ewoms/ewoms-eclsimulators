# .. cmake_module::
#
# This module's content is executed whenever a Dune module requires or
# suggests ewoms-eclsimulators!
#

find_package(Boost REQUIRED)

# either std::any or std::experimental::any needs to be supported
find_package(StdAny REQUIRED)
