cmake_minimum_required(VERSION 2.8)


#######################################
#         PROJECT INFORMATION         #
#-------------------------------------#
# This CMakeLists.txt file is for the #
# CoolProp thermodynamic library      #
# written by Ian Bell. The following  #
# section contains project and        #
# version information.                #
#######################################

# Project name
set(project_name "CoolProp")
set(app_name ${project_name})
project(${project_name})


#######################################
#           BUILD OPTIONS             #
#-------------------------------------#
# These options are available to be   #
# modified in the build process.      #
# packages may want to modify these   #
# to suit, or just leave as defaults. #
#######################################

option (COOLPROP_STATIC_LIBRARY
        "Build and install CoolProp as a STATIC library (.lib, .a) as opposed to SHARED (.dll, .so)"
        ON)

#######################################
#         FIND ALL SOURCES            #
#-------------------------------------#
# The project is organised with       #
# split includes and source folders   #
# this makes it easier for developers #
# to quickly find relevant includes.  #
# This section finds all sources,     #
# headers and corrosponding dirs.     #
#######################################

message(STATUS ${CMAKE_CURRENT_SOURCE_DIR})
file(GLOB APP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/${project_name}/*.cpp")
file(GLOB APP_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/${project_name}/*.h")

set (APP_INCLUDE_DIRS "")
foreach (_headerFile ${APP_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND APP_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES APP_INCLUDE_DIRS)

include_directories(${APP_INCLUDE_DIRS})

#######################################
#         REQUIRED MODULES            #
#-------------------------------------#
# CoolProp requires some standard OS  #
# features, these include:            #
# DL (CMAKE_DL_LIBS) for REFPROP      #
#######################################
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/dev/cmake/Modules/")

if(UNIX)
    find_package (${CMAKE_DL_LIBS} REQUIRED)
endif()

#######################################
#         MAKE ARTIFACTS              #
#-------------------------------------#
# In this section we define the       #
# artifacts (exes, libs) that will be #
# made for coolprop, these include    #
# customisation from earier options.  #
#######################################

###      COOLPROP LIB or DLL        ###
if (COOLPROP_STATIC_LIBRARY)
    add_library(${app_name} STATIC ${APP_SOURCES})
else()
    add_library(${app_name} SHARED ${APP_SOURCES})
endif()
