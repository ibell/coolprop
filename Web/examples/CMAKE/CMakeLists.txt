cmake_minimum_required(VERSION 2.6)

#######################################
#         PROJECT INFORMATION         #
#######################################
set(project_name "MY Project")
set(app_name "my_project")
project(${project_name})


#######################################
#         FIND ALL SOURCES            #
#######################################
# 1) Find all source files
file(GLOB_RECURSE APP_SOURCES "${app_name}/*.cpp")
file(GLOB_RECURSE APP_HEADERS "${app_name}/*.hxx")

# 2) Get all include directories from src
set (APP_INCLUDE_DIRS "")
foreach (_headerFile ${APP_HEADERS})
    get_filename_component(_dir ${_headerFile} PATH)
    list (APPEND APP_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES APP_INCLUDE_DIRS)

# 3) Set executable name and sources
add_executable(${app_name} ${APP_SOURCES})


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/libs/CoolProp/dev/cmake/Modules")

################################################
#        INCLUDED LIBRARIES FROM /LIB          #
################################################
add_subdirectory("libs/CoolProp")
include_directories("libs/CoolProp/CoolProp")
target_link_libraries (${app_name} CoolProp)

if(UNIX)
  find_package (${CMAKE_DL_LIBS} REQUIRED)
  target_link_libraries (${app_name} dl)
endif()



