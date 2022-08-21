set(VORPALINE_ARCH_64 true)
include(${GEOGRAM_SOURCE_DIR}/cmake/platforms/Linux-icc.cmake)

IF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm")
    add_flags(CMAKE_CXX_FLAGS -m64)
    add_flags(CMAKE_C_FLAGS -m64)
ELSE(${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm")
    message("NOT ON ARM CPU! Omitting arm specific flag.")
ENDIF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm")
