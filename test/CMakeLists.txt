#CFLAGS
include_directories(tinytest)
add_executable(test_fdb_internals test.c tinytest/tinytest.c)
target_link_libraries(test_fdb_internals z)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY bin)

add_test(NAME test_fastDBarcode_internals COMMAND test_fdb_internals)
configure_file(integration.py integration.py COPYONLY)
add_test(NAME test_fastDBarcode_output COMMAND python integration.py)
