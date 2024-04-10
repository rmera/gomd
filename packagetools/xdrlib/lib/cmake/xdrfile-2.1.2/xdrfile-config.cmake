# Config file for the INSTALLED package
# Allow other CMake projects to find this package if it is installed
# Requires the use of the standard CMake module CMakePackageConfigHelpers

set ( xdrfile_VERSION 2.1.2 )

set ( xdrfile_FOUND TRUE )

find_path( xdrfile_INCLUDE_DIRS
    NAMES "xdrfile.h" "xdrfile_xtc.h"
    HINTS "/opt/include"
)

find_library( xdrfile_LIBRARIES
    NAMES xdrfile libxdrfile 
    HINTS "/opt/lib" 
)

