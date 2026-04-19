cmake_minimum_required(VERSION 3.16)
project(sde_models CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# =============================================================================
# Matplot++ — скачивается автоматически при первой сборке через FetchContent.
# Не требует ручной установки на целевой машине.
# =============================================================================
include(FetchContent)

FetchContent_Declare(
        matplotplusplus
        GIT_REPOSITORY https://github.com/alandefreitas/matplotplusplus.git
        GIT_TAG        origin/master
        GIT_SHALLOW    TRUE
)

set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BUILD_TESTS    OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(matplotplusplus)

# =============================================================================
# Исполняемый файл
# =============================================================================
add_executable(sde_models src/main.cpp)

target_link_libraries(sde_models PRIVATE matplot)

target_include_directories(sde_models PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Создаём папку plots рядом с бинарником автоматически
add_custom_command(TARGET sde_models POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory
        $<TARGET_FILE_DIR:sde_models>/../plots
        COMMENT "Creating plots/ directory"
)