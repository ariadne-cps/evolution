find_library(ARIADNE_EVOLUTION_LIBRARY NAMES ariadne-evolution)
find_path(ARIADNE_EVOLUTION_INCLUDE_DIR ariadne-evolution.hpp
          PATH_SUFFIXES ariadne-evolution)

find_package(Threads REQUIRED)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AriadneEvolution DEFAULT_MSG
    ARIADNE_EVOLUTION_LIBRARY
    ARIADNE_EVOLUTION_INCLUDE_DIR
)

if(AriadneEvolution_FOUND)
    get_filename_component(ARIADNE_EVOLUTION_INCLUDE_PARENT_DIR
        "${ARIADNE_EVOLUTION_INCLUDE_DIR}" DIRECTORY)

    if(NOT TARGET AriadneEvolution::ariadne-evolution)
        add_library(AriadneEvolution::ariadne-evolution UNKNOWN IMPORTED)
        set_target_properties(AriadneEvolution::ariadne-evolution PROPERTIES
            IMPORTED_LOCATION "${ARIADNE_EVOLUTION_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES
                "${ARIADNE_EVOLUTION_INCLUDE_PARENT_DIR};${ARIADNE_EVOLUTION_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "Threads::Threads"
        )
    endif()
endif()
