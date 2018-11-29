cmake_minimum_required(VERSION 3.6)
project(HyperSolveMiniApp)
set(CMAKE_CXX_STANDARD 11)

set(MINIAPP_HEADERS
        src/Element.h
        src/Ddata.h
        src/HyperSolveMiniApp.h
        src/Point.h
        src/Tensor.h)

add_executable(${PROJECT_NAME} src/HyperSolveMiniApp.cpp ${MINIAPP_HEADERS})

target_include_directories(${PROJECT_NAME} PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        )

install(FILES ${HEADER_FILES}
        DESTINATION include)

install(TARGETS ${PROJECT_NAME}
        DESTINATION bin)
