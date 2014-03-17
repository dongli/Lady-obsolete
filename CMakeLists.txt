# ------------------------------------------------------------------------------
# CMake build script for Lady.
#
# Copyright 2013.
#
# Authors:
#   - Li Dong <dongli@lasg.iap.ac.cn>
# ------------------------------------------------------------------------------

cmake_minimum_required (VERSION 2.8)

# ------------------------------------------------------------------------------
project (lady CXX)

add_definitions (-DDEBUG)
set (CMAKE_CXX_FLAGS "-std=c++11 -Wno-c++11-extensions")

option (USE_SYSTEM_MLPACK "Use system library MLPACK" ON)
option (USE_SYSTEM_NLOPT "Use system library NLopt" ON)

# ------------------------------------------------------------------------------
# external libraries
# MLPACK
set (mlpack_internal_root "${PROJECT_SOURCE_DIR}/external/mlpack-1.0.8")
if (EXISTS "${mlpack_internal_root}" AND NOT ${USE_SYSTEM_MLPACK})
    # use internal MLPACK
    message ("@@ Use internal library MLPACK!")
    add_subdirectory ("${mlpack_internal_root}")
    include_directories ("${PROJECT_BINARY_DIR}/include") # for MLPACK
else ()
    # use user provided MLPACK
    if (DEFINED ENV{MLPACK_ROOT})
        message ("@@ Use user provided library MLPACK.")
        message ("@@ MLPACK_ROOT = $ENV{MLPACK_ROOT}")
        set (MLPACK_INCLUDE_DIR "$ENV{MLPACK_ROOT}/include")
        set (MLPACK_LIBRARY_DIR "$ENV{MLPACK_ROOT}/lib")
        find_package(LibXml2 2.6.0 REQUIRED)
        include_directories (${MLPACK_INCLUDE_DIR} ${LIBXML2_INCLUDE_DIR})
        link_directories (${MLPACK_LIBRARY_DIR})
    else ()
        message (FATAL_ERROR
            "CMake couldn't find library MLPACK! "
            "If it have been installed and you know where it is, "
            "set MLPACK_ROOT (e.g. in .bashrc) to it."
        )
    endif ()
endif ()

# NLopt
set (nlopt_internal_root "${PROJECT_SOURCE_DIR}/external/nlopt-2.4.1")
if (EXISTS "${nlopt_internal_root}" AND NOT ${USE_SYSTEM_NLOPT})
    # use internal NLopt
    message ("@@ Use internal library NLopt!")
    add_subdirectory ("${nlopt_internal_root}")
else ()
    # use user provided NLopt
    if (DEFINED ENV{NLOPT_ROOT})
        message ("@@ Use user provided library NLopt.")
        message ("@@ NLOPT_ROOT = $ENV{NLOPT_ROOT}")
        set (NLOPT_INCLUDE_DIR "$ENV{NLOPT_ROOT}/include")
        set (NLOPT_LIBRARY_DIR "$ENV{NLOPT_ROOT}/lib")
    else ()
        message (FATAL_ERROR
            "CMake couldn't find library NLopt! "
            "If it have been installed and you know where it is, "
            "set NLOPT_ROOT (e.g. in .bashrc) to it."
        )
    endif ()
endif ()

# CGAL
find_package (CGAL)
if (NOT CGAL_FOUND)
    if (DEFINED ENV{CGAL_ROOT})
        message ("@@ Use user provided library CGAL!")
        message ("@@ CGAL_ROOT = $ENV{CGAL_ROOT}")
    else ()
        message (FATAL_ERROR
            "CMake couldn't find library CGAL! "
            "If it have been installed and you know where it is, "
            "set CGAL_ROOT (e.g. in .bashrc) to it."
        )
    endif ()
endif ()
include (${CGAL_USE_FILE})

# Boost
find_package (Boost REQUIRED)
if (NOT Boost_FOUND)
    if (DEFINED ENV{Boost_ROOT})
        message ("@@ Use user provided library Boost!")
        message ("@@ Boost_ROOT = $ENV{Boost_ROOT}")
    else ()
        message (FATAL_ERROR
            "CMake couldn't find library Boost! "
            "If it have been installed and you know where it is, "
            "set Boost_ROOT (e.g. in .bashrc) to it."
        )
    endif ()
endif ()

# GeoMTK
add_subdirectory ("${PROJECT_SOURCE_DIR}/external/geomtk")
include_directories (
    ${NLOPT_INCLUDE_DIR}
    ${GEOMTK_INCLUDE_DIRS}
)
link_directories (
    ${NLOPT_LIBRARY_DIR}
    ${GEOMTK_LIBRARY_DIRS}
)

# ------------------------------------------------------------------------------
# source directory structure
set (source_directories
    "${PROJECT_SOURCE_DIR}/src"
    "${PROJECT_SOURCE_DIR}/src/Utils"
    "${PROJECT_SOURCE_DIR}/src/Parcel"
    "${PROJECT_SOURCE_DIR}/src/Advection"
    "${PROJECT_SOURCE_DIR}/src/Advection/Tracer"
    "${PROJECT_SOURCE_DIR}/src/Advection/TestCase"
    "${PROJECT_SOURCE_DIR}/src/Advection/ModelAdaptor"
)

# ------------------------------------------------------------------------------
# collect sources and headers
foreach (dir ${source_directories})
    include_directories ("${dir}")
    # header files
    file (GLOB tmp1 "${dir}/*.h")
    list (APPEND headers ${tmp1})
    # source files
    aux_source_directory ("${dir}" tmp2)
    list (APPEND sources ${tmp2})
endforeach ()

# ------------------------------------------------------------------------------
# library targets
add_library (lady STATIC ${headers} ${sources})
target_link_libraries (lady mlpack nlopt CGAL)
add_dependencies (lady geomtk)

# ------------------------------------------------------------------------------
# executable targets
foreach (dir ${source_directories})
    if (EXISTS "${dir}/demo")
        include_directories ("${dir}/demo")
        aux_source_directory ("${dir}/demo" demo_paths)
        foreach (demo_path ${demo_paths})
            get_filename_component (demo ${demo_path} NAME_WE)
            add_executable (${demo} ${demo_path})
            target_link_libraries (${demo} geomtk lady)
        endforeach ()
    endif ()
endforeach ()

# ------------------------------------------------------------------------------
# install rules
if (${CMAKE_GENERATOR} MATCHES "Unix Makefiles")
    install (TARGETS lady ARCHIVE DESTINATION lib)
    foreach (header ${headers})
        if (${header} MATCHES "lady.h")
            install (FILES ${header}
                DESTINATION "include"
            )
        else ()
            install (FILES ${header}
                DESTINATION "include/lady"
            )
        endif ()
    endforeach ()
endif ()
