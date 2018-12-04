cmake_minimum_required(VERSION 3.6)
project(HyperSolveMiniApp)
set(CMAKE_CXX_STANDARD 11)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -Wno-int-in-bool-context")
endif()

set(MINIAPP_HEADERS
        src/Element.h
        src/Ddata.h
        src/HyperSolveMiniApp.h
        src/Point.h
        src/TetMesh.h
        src/EdgeGradients.h)

find_package(Eigen3 REQUIRED)

add_executable(${PROJECT_NAME} src/HyperSolveMiniApp.cpp ${MINIAPP_HEADERS})

target_include_directories(${PROJECT_NAME} PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${EIGEN3_INCLUDE_DIRS}
        )

install(FILES ${HEADER_FILES}
        DESTINATION include)

install(TARGETS ${PROJECT_NAME}
        DESTINATION bin)
