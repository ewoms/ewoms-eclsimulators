# .. cmake_module::
#
# This module's content is executed whenever a Dune module requires or
# suggests ewoms-eclsimulators!
#

find_package(Boost
  COMPONENTS filesystem regex system date_time
  REQUIRED)