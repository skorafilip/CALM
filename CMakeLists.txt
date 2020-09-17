#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required (VERSION 2.6 FATAL_ERROR)
project (CALM)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
find_package(ROOT REQUIRED)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include/
                    ${ROOT_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cxx)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(calm ${sources} ${headers})
target_link_libraries(calm ${ROOT_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory
#
set(CALM_SCRIPTS
    share/particles.data
    events.ini
    config.ini
    run_CALM.sh
    Doxyfile)
foreach(_script ${CALM_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Makeing events directory
#
make_directory(events)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS calm DESTINATION bin)


