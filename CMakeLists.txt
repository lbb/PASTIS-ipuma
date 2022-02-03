cmake_minimum_required(VERSION 3.11)
project(pastis LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/Modules)

include_directories("plog/include")
include_directories(
        ${PROJECT_SOURCE_DIR}/seqan/include 
        ${PROJECT_SOURCE_DIR}/CombBLAS/install/include
        ${PROJECT_SOURCE_DIR}/CombBLAS/Applications
        ${PROJECT_SOURCE_DIR}/src/
        ${PROJECT_SOURCE_DIR}/src/ipuma-sw/
        )

link_directories(
    ${PROJECT_SOURCE_DIR}/CombBLAS/install/lib
    )

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DSEQAN_ENABLE_DEBUG=0")
add_compile_definitions(NDEBUG)

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

FILE(GLOB PASTIS_CXX_SRCS src/*.cpp src/kmer/*.cpp src/align/*.cpp src/ipuma-sw/*.cpp)
FILE(GLOB PASTIS_CXX_INCS inc/*.hpp inc/kmer/*.hpp inc/align/*.hpp src/ipuma-sw/*.hpp src/ipuma-sw/*.h)
add_executable(pastis
        ${PASTIS_CXX_SRCS}
        ${PASTIS_CXX_INCS}
        )

# target_link_libraries(pastis PRIVATE -fsanitize=address)
target_link_libraries(pastis PUBLIC CombBLAS GraphGenlib Usortlib)
if (MPI_FOUND)
        message(STATUS "MPI found.")
    include_directories(${MPI_INCLUDE_PATH})
    target_link_libraries(pastis PUBLIC ${MPI_LIBRARIES})
endif()


add_definitions(-DPOPLAR_ENABLED)
include(AddGPLibrary)
message(STATUS "Building for IPU with POPLAR")
set(POPLIBS_ENABLED_IPU_ARCH_NAMES "cpu;ipu2")
# if(NOT DEFINED POPLIBS_ENABLED_IPU_ARCH_NAMES)
#     set(POPLIBS_ENABLED_IPU_ARCH_NAMES "cpu;ipu1;ipu2" CACHE STRING
#         "The architectures to build the poplibs library objects for")
# endif()
message("Enabled popc architectures: ${POPLIBS_ENABLED_IPU_ARCH_NAMES}")

set(DEFAULT_TEST_VARIANTS "cpu;IpuModel2;")

set(POPC_EXECUTABLE "popc")

# Set up popc flags.
set(POPC_FLAGS -I ${CMAKE_SOURCE_DIR}/include)

if(CMAKE_BUILD_TYPE STREQUAL Debug)
  list(APPEND POPC_FLAGS -g)
elseif(CMAKE_BUILD_TYPE STREQUAL RelWithDebInfo)
  list(APPEND POPC_FLAGS -O3 -g)
else()
  list(APPEND POPC_FLAGS -O3)
endif()

# If requested, optimise the debug builds too.
if (OPTIMISE_IPU_DEBUG AND CMAKE_BUILD_TYPE STREQUAL Debug)
  list(APPEND POPC_FLAGS -O3)
endif()

add_subdirectory(src/ipuma-sw)
# add_library(IPUMA_LIBRARY "ipuma-sw/matrix.cpp" "ipuma-sw/similarity.cpp" "ipuma-sw/encoding.cpp")

# include logging library

set(IPUMA_LINK_LIBRARIES ${IPUMA_LINK_LIBRARIES} IPUMA_LIBRARY poplar)
target_link_libraries(pastis PUBLIC ${IPUMA_LINK_LIBRARIES})


if(OPENMP_FOUND OR OpenMP_CXX_FOUND)
    # Set THREADED if OpenMP is found
    message("Threaded")
    target_compile_definitions(pastis PUBLIC THREADED)
    if(TARGET OpenMP::OpenMP_CXX)
        target_link_libraries(pastis PUBLIC OpenMP::OpenMP_CXX)
    else()
        target_compile_options(pastis PUBLIC "${OpenMP_CXX_FLAGS}")
        target_link_libraries(pastis PUBLIC "${OpenMP_CXX_FLAGS}")
    endif()
endif()
