# ------------------------------------------------------------------------------
# CMake build script for Lady.
#
# Copyright 2013.
#
# Authors:
#   - Li Dong <dongli@lasg.iap.ac.cn>
# ------------------------------------------------------------------------------

cmake_minimum_required (VERSION 2.8)

project (lady CXX)

# ------------------------------------------------------------------------------
# external libraries
find_package (Armadillo)
if (NOT ARMADILLO_FOUND)
    if (DEFINED ENV{ARMADILLO_ROOT})
        message ("@@ Use user provided library Armadillo!")
        message ("@@ ARMADILLO_ROOT = $ENV{ARMADILLO_ROOT}")
        find_package (Armadillo HINTS $ENV{ARMADILLO_ROOT})
    else ()
        message (FATAL_ERROR
            "CMake couldn't find library Armadillo! "
            "If it have been installed and you know where it is, "
            "set ARMADILLO_ROOT (e.g. in .bashrc) to it."
        )
    endif ()
endif ()

add_subdirectory (${PROJECT_SOURCE_DIR}/external/geomtk)

# ------------------------------------------------------------------------------
# source directory structure
set (source_directories
    "${PROJECT_SOURCE_DIR}/src"
    "${PROJECT_SOURCE_DIR}/src/Advection"
    "${PROJECT_SOURCE_DIR}/src/Advection/Tracer"
    "${PROJECT_SOURCE_DIR}/src/Advection/TestCase"
)

# ------------------------------------------------------------------------------
# collect sources and headers
foreach (dir ${source_directories})
    include_directories ("${dir}")
    # source files
    aux_source_directory ("${dir}" tmp1)
    list (APPEND sources ${tmp1})
    # header files
    file (GLOB tmp2 "${dir}/*.h")
    list (APPEND headers ${tmp2})
endforeach ()

# ------------------------------------------------------------------------------
# library targets
add_library (lady STATIC ${sources})
include_directories (
    ${ARMADILLO_INCLUDE_DIRS}
    ${GEOMTK_INCLUDE_DIRS}
)
link_directories (
    ${GEOMTK_LIBRARY_DIRS}
)
target_link_libraries (lady armadillo)
add_dependencies (lady geomtk)

install (TARGETS lady
    ARCHIVE DESTINATION lib
)
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

# ------------------------------------------------------------------------------
# executable targets
foreach (dir ${source_directories})
    if (EXISTS "${dir}/demo")
        include_directories ("${dir}/demo")
        aux_source_directory ("${dir}/demo" demo_paths)
        foreach (demo_path ${demo_paths})
            get_filename_component (demo ${demo_path} NAME_WE)
            add_executable (${demo} ${demo_path})
            target_link_libraries (${demo} armadillo geomtk lady)
        endforeach ()
    endif ()
endforeach ()

# ------------------------------------------------------------------------------
# install rules
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
