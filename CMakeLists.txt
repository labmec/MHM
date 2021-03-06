# @file MHM/CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 2.8.0)

set(CMAKE_CONFIGURATION_TYPES Debug Release CACHE STRING  INTERNAL FORCE )

project(MHM)

IF (WIN32)
	# Prevents timespec redefinition problem with Visual Studio 2015
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_STRUCT_TIMESPEC")
ELSE (WIN32)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")#this flag is for preventing errors when casting from double to float
	if ( NOT APPLE AND NOT CMAKE_BUILD_TYPE )
		set (CMAKE_BUILD_TYPE Release CACHE STRING "Options: Release Debug None MinSizeRel RelWithDebInfo" FORCE) 
	endif()
ENDIF()


set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

IF (WIN32)
	find_package( PZ PATHS "C:/Arquivos de Programas/PZ" REQUIRED NO_DEFAULT_PATH)
ELSE()
	find_package( PZ REQUIRED )
ENDIF (WIN32)


IF(APPLE)
   # Mac OS X specific code
   add_definitions(-DMACOSX)
   set(MACOSX MACOSX)
#   SET(OperatingSystem "Mac OS X")
ENDIF(APPLE)

include_directories(${PZ_INCLUDE_DIRS})
include_directories(${PZ_PTHREAD_INCLUDE})

# Location to search files or libraries
set(SEARCH_DIRS ENV PATH /usr/local/opt/ /usr/local/opt/qt5/lib/cmake/ /usr/include /usr/local/include /usr/local/include/boost /usr/local/include/log4cxx /usr/local/include/tbb /usr/lib /usr/lib64 /usr/local/lib ./. ./.. ./../libs ./../libs/fad ./../libs/metis ./../libs/metis/include ./../libs/metis/lib ./../externallibs ./../externallibs/fad ./../externallibs/pthread ./../externallibs/metis ./../externallibs/metis/include ./../externallibs/metis/lib ./..metis ./../metis/include ./../metis/lib ./externallibs/lib ./externallibs/include ./../externallibs/lib ./../externallibs/include)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	set(PZDEBUG PZDEBUG)
	add_definitions(-DPZDEBUG)
endif()
if (CMAKE_BUILD_TYPE STREQUAL "Release")
	set(NODEBUG NODEBUG)
	set(ZERO_INTERNAL_RESIDU ZERO_INTERNAL_RESIDU)
    set(MAKEINTERNAL MAKEINTERNAL)
	add_definitions(-DNODEBUG)
endif()


if(WIN32)
	#add a define saying it is a VC compiler
	set(VC "#Will add a define VC on pz_config.h")
	add_definitions(-DVC)
	#define _USE_MATH_DEFINES for example: M_PI constant definitions.
	add_definitions(-D_USE_MATH_DEFINES)
	#disabling VC warnings
	add_definitions(-D_SCL_SECURE_NO_WARNINGS)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    #enabling /bigobj
    add_definitions("/bigobj")

    FOREACH(FLAG_TYPE EXE MODULE SHARED)
        STRING (REPLACE "INCREMENTAL:YES" "INCREMENTAL" FLAG_TMP "${CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG}")
        STRING (REPLACE "INCREMENTAL:NO" "INCREMENTAL" FLAG_TMP ${FLAG_TMP})
        STRING (REPLACE "INCREMENTAL" "INCREMENTAL:NO" FLAG_TMP ${FLAG_TMP})
        STRING (REPLACE "/EDITANDCONTINUE" "" FLAG_TMP ${FLAG_TMP})
        SET(CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG "${FLAG_TMP}" CACHE STRING "Overriding default debug ${FLAG_TYPE} linker flags." FORCE)
        MARK_AS_ADVANCED (CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG)
    ENDFOREACH ()

endif()


# Boost library dependencies
IF (WIN32)
	set(Boost_USE_STATIC_LIBS   ON)
	set(Boost_USE_MULTITHREADED ON)
        set(CMAKE_PREFIX_PATH "C:/local/boost_1_72_0" ${CMAKE_PREFIX_PATH})
ENDIF()

find_package( Boost 1.71 COMPONENTS unit_test_framework date_time atomic REQUIRED )


add_executable(MHM main.cpp TPZMatLaplacianHybrid.h TPZMatLaplacianHybrid.cpp meshgen.h meshgen.cpp)

target_link_libraries(MHM pz )

add_executable(MHM_HDiv main_MHM_HDiv.cpp TPZMatLaplacianHybrid.h TPZMatLaplacianHybrid.cpp meshgen.h meshgen.cpp)

target_link_libraries(MHM_HDiv pz )

add_executable(MHM_HDiv_elast main_MHM_HDiv_elast.cpp TPZMatLaplacianHybrid.h TPZMatLaplacianHybrid.cpp meshgen.h meshgen.cpp)

target_link_libraries(MHM_HDiv_elast pz )

#add_executable(MHMMaze mainMaze.cpp TPZMatLaplacianHybrid.h TPZMatLaplacianHybrid.cpp meshgen.h meshgen.cpp ConfigCasesMaze.h)

#target_link_libraries(MHMMaze pz)

add_executable(MHM_HDiv_elastHetero main_MHM_HDiv_elastHetero.cpp TPZMatLaplacianHybrid.h TPZMatLaplacianHybrid.cpp meshgen.h meshgen.cpp)

target_link_libraries(MHM_HDiv_elastHetero pz )

add_executable(MHMElast mainElast.cpp TPZMatLaplacianHybrid.h TPZMatLaplacianHybrid.cpp meshgen.h meshgen.cpp)

target_link_libraries(MHMElast pz)

add_executable(ElasticidadeHetero heteroElast.cpp TPZMatLaplacianHybrid.h TPZMatLaplacianHybrid.cpp meshgen.h meshgen.cpp)

target_link_libraries(ElasticidadeHetero pz)


configure_file(young_modulus.txt young_modulus.txt COPYONLY)
configure_file(Poisson_coef.txt Poisson_coef.txt COPYONLY)


#add_executable(MHM_Compare MHM_Compare.cpp TPZMatLaplacianHybrid.h TPZMatLaplacianHybrid.cpp meshgen.h meshgen.cpp)

#target_link_libraries(MHM_Compare pz)

file(COPY mhm_hdiv_heterogeneous_elast_exec.txt DESTINATION .)
file(COPY mhmelast_exec.txt DESTINATION .)
file(COPY Data_13_Set DESTINATION .)

#Enabling nightly tests
ENABLE_TESTING()
include(CTest)
