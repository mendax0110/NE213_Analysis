cmake_minimum_required(VERSION 3.16)
project(NE213_Analysis VERSION 1.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Qt6 REQUIRED COMPONENTS Core Widgets Charts)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic -fPIC)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    add_compile_options(/W4)
endif()

set(CMAKE_AUTOMOC ON)

set(SOURCES
        main.cpp
        src/WaveformData.cpp
        src/WaveformProcessor.cpp
        src/PSDAnalyzer.cpp
        src/Statistics.cpp
        src/QtPlotter.cpp
)

set(HEADERS
        include/WaveformData.h
        include/WaveformProcessor.h
        include/PSDAnalyzer.h
        include/Statistics.h
        include/QtPlotter.h
)

add_executable(${PROJECT_NAME}
        ${SOURCES}
        ${HEADERS}
)

target_include_directories(${PROJECT_NAME} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_BINARY_DIR}
)

target_link_libraries(${PROJECT_NAME} PRIVATE
        Qt6::Core
        Qt6::Widgets
        Qt6::Charts
)

set_target_properties(${PROJECT_NAME} PROPERTIES
        AUTOMOC ON
        AUTOUIC ON
        AUTORCC ON
)

install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION bin
)