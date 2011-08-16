# - Try to find Eigen2 lib
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(Eigen2 2.0.3)
# to require version 2.0.3 to newer of Eigen2.
#
# Once done this will define
#
#  EIGEN2_FOUND - system has eigen lib with correct version
#  EIGEN2_INCLUDE_DIR - the eigen include directory
#  EIGEN2_VERSION - eigen version

# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
# Redistribution and use is allowed according to the terms of the BSD license.

if (EIGEN2_INCLUDE_DIR)

  # in cache already
  set(EIGEN2_FOUND TRUE)

else (EIGEN2_INCLUDE_DIR)

find_path(EIGEN2_INCLUDE_DIR NAMES Eigen/Core
     PATHS
     ${INCLUDE_INSTALL_DIR}
     ${KDE4_INCLUDE_DIR}
     PATH_SUFFIXES eigen2
   )

if(EIGEN2_INCLUDE_DIR)
  set(EIGEN2_FOUND TRUE)
endif(EIGEN2_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Eigen2 DEFAULT_MSG EIGEN2_INCLUDE_DIR EIGEN2_VERSION_OK)

mark_as_advanced(EIGEN2_INCLUDE_DIR)

endif(EIGEN2_INCLUDE_DIR)

