#############################################################
# CMake requirement (From Modern CMake)
#############################################################
cmake_minimum_required(VERSION 3.20)

get_filename_component( TARGET_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME )
message( STATUS "\n-------------  Config ${PACKAGE_NAME}::${TARGET_NAME}: ${CMAKE_CURRENT_SOURCE_DIR}  -------------" )

file(GLOB SOURCES src/*.cpp     src/*.c     )
file(GLOB HEADERS src/*.h       src/*.hpp   )
file(GLOB INTERFACE_HEADERS include/${TARGET_NAME}/*.h include/${TARGET_NAME}/*.hpp )

cmake_path(SET abs_include_dir NORMALIZE "${CMAKE_CURRENT_SOURCE_DIR}/include")

if(NOT PCL_DIR)
    get_filename_component(_smrobot_default_pcl_dir
        "${CMAKE_SOURCE_DIR}/../ThirdPartyLib/PCL1.12.1/install/cmake"
        ABSOLUTE)
    if(EXISTS "${_smrobot_default_pcl_dir}/PCLConfig.cmake")
        set(PCL_DIR "${_smrobot_default_pcl_dir}" CACHE PATH
            "Directory containing PCLConfig.cmake")
    endif()
endif()
find_package(PCL 1.12 REQUIRED COMPONENTS
    common filters kdtree search features surface)
get_filename_component(_smrobot_pcl_install_root "${PCL_DIR}" DIRECTORY)

add_library( ${TARGET_NAME} STATIC ${SOURCES} ${HEADERS} ${INTERFACE_HEADERS} )
add_library( ${PACKAGE_NAME}::${TARGET_NAME} ALIAS ${TARGET_NAME})
target_compile_definitions(${TARGET_NAME} PRIVATE HAVE_QHULL)
if(MSVC)
    target_compile_options(${TARGET_NAME} PRIVATE /UDISABLE_QHULL)
endif()

target_include_directories(${TARGET_NAME}
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${PACKAGE_NAME}/include>
    PRIVATE
        $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
        "${_smrobot_pcl_install_root}/3rdParty/FLANN/include"
        "${_smrobot_pcl_install_root}/3rdParty/Qhull/include"
)

target_link_libraries ( ${TARGET_NAME}
    PUBLIC ${${PACKAGE_NAME}${TARGET_NAME}_RequiredLibsPublic}
    PRIVATE
        ${${PACKAGE_NAME}${TARGET_NAME}_RequiredLibsPrivate}
        pcl_common
        pcl_filters
        pcl_kdtree
        pcl_search
        pcl_features
        pcl_surface
)

set_target_properties( ${TARGET_NAME} PROPERTIES
    FOLDER ${PACKAGE_NAME}
    DEBUG_POSTFIX _shared_${CMAKE_DEBUG_POSTFIX}
    RELEASE_POSTFIX _shared_${CMAKE_RELEASE_POSTFIX}
)

if(BuildTests)
    add_subdirectory(regression/AlgorithmReproductionTest)
    add_subdirectory(regression/DemoThicknessPredictionTest)
    add_subdirectory(regression/RotationalSurfaceFitterTest)
endif()

if( ExportPacakge${PACKAGE_NAME}_${TARGET_NAME} )
    function_InstallTarget(
        TARGET_NAME ${TARGET_NAME}
        PACKAGE_NAME ${PACKAGE_NAME}
        PACKAGE_INTERFACE_HEADERS ${abs_include_dir}
    )
endif()

message( STATUS "-------------  Config ${PACKAGE_NAME}::${TARGET_NAME}: ${CMAKE_CURRENT_SOURCE_DIR}  -------------\n" )
