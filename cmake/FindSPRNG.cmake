#
# this module look for SPRNG (http://hdf.ncsa.uiuc.edu) support
# it will define the following values
#
# SPRNG_INCLUDE_PATH = where sprng.h can be found
# SPRNG_LIBRARY = the library to link against (sprng etc)
#

SET(TRIAL_LIBRARY_PATHS
  /usr/local/lib
  /sw/lib
  ${CMAKE_SOURCE_DIR}/lib
  $ENV{SPRNG_HOME}/lib
)

SET(TRIAL_INCLUDE_PATHS
  /usr/include/sprng/
  /usr/local/include
  /sw/include
  ${CMAKE_SOURCE_DIR}/include
  $ENV{SPRNG_HOME}/include
)

find_library(SPRNG_LIBRARY sprng ${TRIAL_LIBRARY_PATHS})
find_path(SPRNG_INCLUDE_DIR sprng.h ${TRIAL_INCLUDE_PATHS})
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SPRNG DEFAULT_MSG SPRNG_LIBRARY SPRNG_INCLUDE_DIR)

mark_as_advanced(
SPRNG_LIBRARY
SPRNG_INCLUDE_DIR
)
