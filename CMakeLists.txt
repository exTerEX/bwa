cmake_minimum_required(VERSION 3.20)
project(bwa LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

enable_testing()

add_library(bwa_core
    src/fasta.cpp
    src/index_io.cpp
)

target_include_directories(bwa_core PUBLIC include)

add_executable(bwa
    src/main.cpp
)
target_link_libraries(bwa PRIVATE bwa_core)

add_executable(bwa_tests
    tests/test_fm_index.cpp
)
target_link_libraries(bwa_tests PRIVATE bwa_core)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(bwa_core PRIVATE
        -O3
        -Wall
        -Wextra
        -Wpedantic
        -Wconversion
    )
    target_compile_options(bwa PRIVATE
        -O3
        -Wall
        -Wextra
        -Wpedantic
        -Wconversion
    )
    target_compile_options(bwa_tests PRIVATE
        -O3
        -Wall
        -Wextra
        -Wpedantic
        -Wconversion
    )
endif()

add_test(NAME bwa_tests COMMAND bwa_tests)
