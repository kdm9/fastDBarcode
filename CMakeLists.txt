cmake_minimum_required(VERSION 2.6)
project(fastDBarcode)

# Cmake options
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
enable_testing()

# CC options
add_definitions(-std=gnu99 -O3 -fgnu89-inline)
include_directories(${fastDBarcode_SOURCE_DIR}/src)
link_directories(${fastDBarcode_BINARY_DIR}/lib)

# Targets
add_executable(fastDBarcode src/main.c)
target_link_libraries(fastDBarcode z)


INSTALL(TARGETS fastDBarcode DESTINATION "bin")
