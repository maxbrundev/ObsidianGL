project(ObsidianGL LANGUAGES CXX)

file(GLOB_RECURSE OBSIDIAN_FILES CONFIGURE_DEPENDS
    "src/*.cpp"
    "include/*.h"
)

add_library(ObsidianGL SHARED ${OBSIDIAN_FILES})

target_compile_definitions(ObsidianGL PRIVATE OBSIDIANGL_BUILD)

if(OBSIDIAN_USE_AVX2)
    if(MSVC)
        target_compile_options(ObsidianGL PRIVATE /arch:AVX2)
        target_compile_definitions(ObsidianGL PRIVATE __AVX2__)
    else()
        target_compile_options(ObsidianGL PRIVATE -mavx2)
    endif()
endif()

target_include_directories(ObsidianGL
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_link_libraries(ObsidianGL 
    PUBLIC 
        glm
)

set_target_properties(ObsidianGL PROPERTIES
    OUTPUT_NAME "ObsidianGL"
    FOLDER "Libraries"
    OBJECT_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/Bin-Int/ObsidianGL
)