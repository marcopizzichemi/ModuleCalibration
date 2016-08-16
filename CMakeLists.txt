#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(ModuleCalibration)

# Load some basic macros which are needed later on
include(FindROOT.cmake)

#find_package(OpenMP)

#if(OPENMP_FOUND)
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

# Create a library called "libEvent" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

set(INCLUDE_DIRECTORIES
${ROOT_INCLUDE_DIR} 
)

include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
)

#set (CMAKE_CXX_FLAGS "-Wshadow -Wvla -Wsign-compare")


#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
# option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
# if(WITH_GEANT4_UIVIS)
#   find_package(Geant4 REQUIRED ui_all vis_all)
# else()
#   find_package(Geant4 REQUIRED)
# endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
# include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include )
#                     ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
#ROOT_LIBRARIES := $(filter-out -lNew,$(ROOT_LIBRARIES))

add_executable(ModuleCalibration ModuleCalibration.cc ${sources} ${headers})
target_link_libraries(ModuleCalibration ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} -lSpectrum -lMLP -lTreePlayer)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build OpNovice. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(g4matrix_SCRIPTS
#    g4matrix.out
#    g4matrix.in
#    optPhoton.mac
#     gui.mac
#     icons.mac
#     run.png
#     vis.mac
#     run.mac
#     runSingle.mac
#     matrix.cfg
  )

foreach(_script ${g4matrix_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ModuleCalibration DESTINATION bin)
set(CMAKE_BUILD_TYPE Release)
# set (CMAKE_CXX_FLAGS "-Wvla -Wsign-compare")
set (CMAKE_CXX_FLAGS "-W -Wall -pedantic -Wno-non-virtual-dtor -Wno-long-long -Wwrite-strings -Wpointer-arith -Woverloaded-virtual -Wno-variadic-macros -pipe -std=c++98")
