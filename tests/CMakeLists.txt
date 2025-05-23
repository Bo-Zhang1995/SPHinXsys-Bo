option(SPHINXSYS_BUILD_2D_EXAMPLES "SPHINXSYS_BUILD_2D_EXAMPLES" ON)
option(SPHINXSYS_BUILD_3D_EXAMPLES "SPHINXSYS_BUILD_3D_EXAMPLES" ON)
option(SPHINXSYS_BUILD_OPTIMIZATION_EXAMPLES "SPHINXSYS_BUILD_OPTIMIZATION_EXAMPLES" ON)
option(SPHINXSYS_BUILD_UNIT_TESTS "SPHINXSYS_BUILD_UNIT_TESTS" ON)
option(SPHINXSYS_BUILD_USER_EXAMPLES "SPHINXSYS_BUILD_USER_EXAMPLES" ON)
option(SPHINXSYS_BUILD_MODULES "SPHINXSYS_BUILD_MODULES" ON)
option(SPHINXSYS_BUILD_PYTHON_INTERFACE "SPHINXSYS_BUILD_PYTHON_INTERFACE" ON)

find_package(GTest CONFIG REQUIRED)
include(GoogleTest)

target_link_libraries(sphinxsys_core INTERFACE GTest::gtest GTest::gtest_main)

if(SPHINXSYS_2D AND SPHINXSYS_BUILD_2D_EXAMPLES)
    ADD_SUBDIRECTORY(2d_examples)
endif()

if(SPHINXSYS_BUILD_OPTIMIZATION_EXAMPLES)
    ADD_SUBDIRECTORY(optimization)
endif()

if(SPHINXSYS_BUILD_PYTHON_INTERFACE)
    ADD_SUBDIRECTORY(test_python_interface)
endif()

if(SPHINXSYS_BUILD_USER_EXAMPLES)
    ADD_SUBDIRECTORY(extra_source_and_tests)
endif()

if(SPHINXSYS_BUILD_UNIT_TESTS)
    ADD_SUBDIRECTORY(unit_tests_src)
endif()

if(SPHINXSYS_USE_SYCL)
    ADD_SUBDIRECTORY(tests_sycl)
endif()

if(SPHINXSYS_BUILD_MODULES)
    add_subdirectory(modules)
endif()

if(SPHINXSYS_3D AND SPHINXSYS_BUILD_3D_EXAMPLES)
    ADD_SUBDIRECTORY(3d_examples)

    if(EMSCRIPTEN)
        ADD_SUBDIRECTORY(webassembly_models)
    endif()
endif()
