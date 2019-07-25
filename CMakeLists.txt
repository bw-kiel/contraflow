cmake_minimum_required(VERSION 2.6)
project(Contraflow)

option(GTEST "Use google test" ON)

include_directories(
	"${PROJECT_BINARY_DIR}"
	"${PROJECT_SOURCE_DIR}/src"
	"${PROJECT_SOURCE_DIR}/src/matrix_stru3"
)


set(CMAKE_CXX_FLAGS "-Wall")
include(CheckCXXCompilerFlag)
	CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
	CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)

if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()


add_subdirectory(src)

enable_testing()


if(GTEST)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lpthread")
        add_subdirectory("${PROJECT_SOURCE_DIR}/ext/gtest-1.8.0")

        enable_testing()
        find_package(GTest REQUIRED)

        include_directories(
                "${gtest_SOURCE_DIR}"
                "${gtest_SOURCE_DIR}/include"
        )

        add_executable(run_tests gtest/allTests.cpp)
        target_link_libraries(run_tests gtest gtest_main contraflow)


endif(GTEST)
