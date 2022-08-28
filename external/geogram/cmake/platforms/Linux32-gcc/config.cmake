set(VORPALINE_ARCH_32 true)
include(${GEOGRAM_SOURCE_DIR}/cmake/platforms/Linux-gcc.cmake)

# m64/32 are i386 specific flags that shouldn't be used on arm builds.
# Deactivating NATIVE_SSE is an optional switch to remove them.
if(NOT ${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm" AND ${NATIVE_SSE})
    add_flags(CMAKE_CXX_FLAGS -m32)
    add_flags(CMAKE_C_FLAGS -m32)
endif()

# Configure FPU to use SSE instructions (IEEE rounding semantics)
# In the default 387 mode, rounding is unpredictable
add_definitions(-mfpmath=sse)

