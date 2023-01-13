set(VORPALINE_ARCH_64 true)
set(VORPALINE_BUILD_DYNAMIC true)
include(${GEOGRAM_SOURCE_DIR}/cmake/platforms/Linux-gcc5.cmake)

# m64/32 are i386 specific flags that shouldn't be used on arm builds.
# Deactivating NATIVE_SSE is an optional switch to remove them.
if(NOT ${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm" AND ${NATIVE_SSE})
    add_flags(CMAKE_CXX_FLAGS -m64)
    add_flags(CMAKE_C_FLAGS -m64)
endif()
