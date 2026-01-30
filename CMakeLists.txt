cmake_minimum_required(VERSION 3.20)

project(Obsidian LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(CMAKE_SUPPRESS_REGENERATION true)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "CMake")

set(SDK_DIR ${CMAKE_SOURCE_DIR}/Build)
set(INT_DIR ${CMAKE_SOURCE_DIR}/Bin-Int)

foreach(CONFIG Debug Release RelWithDebInfo MinSizeRel)
    string(TOUPPER ${CONFIG} CFG)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${CFG} ${SDK_DIR}/bin/${CONFIG})
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${CFG} ${SDK_DIR}/bin/${CONFIG})
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${CFG} ${SDK_DIR}/lib/${CONFIG})
endforeach()

set(SDL_SHARED OFF CACHE BOOL "" FORCE)
set(SDL_STATIC ON  CACHE BOOL "" FORCE)
set(SDL_TEST   OFF CACHE BOOL "" FORCE)
set(SDL_TEST_LIBRARY OFF CACHE BOOL "" FORCE)

add_subdirectory(Dependencies/SDL3)

get_property(SDL_TARGETS DIRECTORY Dependencies/SDL3 PROPERTY BUILDSYSTEM_TARGETS)
foreach(TARGET_NAME IN LISTS SDL_TARGETS)
    set_target_properties(${TARGET_NAME} PROPERTIES FOLDER "Dependencies/SDL3")
endforeach()

add_library(glm INTERFACE)
target_include_directories(glm INTERFACE ${CMAKE_SOURCE_DIR}/Dependencies/glm/include)

add_library(stb_image INTERFACE)
target_include_directories(stb_image INTERFACE ${CMAKE_SOURCE_DIR}/Dependencies/stb_image/include)

add_subdirectory(Sources/ObsidianGL ${CMAKE_SOURCE_DIR}/Sources/ObsidianGL)
add_subdirectory(Sources/OnyxEngine ${CMAKE_SOURCE_DIR}/Sources/OnyxEngine)

set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT OnyxEngine)