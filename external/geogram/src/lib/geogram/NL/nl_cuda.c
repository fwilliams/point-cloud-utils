/*
 *  Copyright (c) 2004-2017, Bruno Levy
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: Bruno Levy
 *
 *     levy@loria.fr
 *
 *     ALICE Project
 *     LORIA, INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 */

#include "nl_cuda.h"
#include "nl_context.h"

/**
 * \file nl_cuda.c
 * \brief Weak-coupling adapter to call CUDA from OpenNL.
 */

/**********************************************************/
/*      CUDA structures and functions                     */
/* Repeated here so that one can compile OpenNL without   */
/* requiring CUDA to be installed in the system.          */
/**********************************************************/

struct cudaDeviceProp {
    char name[256];
    size_t totalGlobalMem;
    size_t sharedMemPerBlock;
    int regsPerBlock;
    int warpSize;
    size_t memPitch;
    int maxThreadsPerBlock;
    int maxThreadsDim[3];
    int maxGridSize[3];
    int clockRate;
    size_t totalConstMem;
    int major;
    int minor;
    size_t textureAlignment;
    size_t texturePitchAlignment;
    int deviceOverlap;
    int multiProcessorCount;
    int kernelExecTimeoutEnabled;
    int integrated;
    int canMapHostMemory;
    int computeMode;
    int maxTexture1D;
    int maxTexture1DMipmap;
    int maxTexture1DLinear;
    int maxTexture2D[2];
    int maxTexture2DMipmap[2];
    int maxTexture2DLinear[3];
    int maxTexture2DGather[2];
    int maxTexture3D[3];
    int maxTexture3DAlt[3];
    int maxTextureCubemap;
    int maxTexture1DLayered[2];
    int maxTexture2DLayered[3];
    int maxTextureCubemapLayered[2];
    int maxSurface1D;
    int maxSurface2D[2];
    int maxSurface3D[3];
    int maxSurface1DLayered[2];
    int maxSurface2DLayered[3];
    int maxSurfaceCubemap;
    int maxSurfaceCubemapLayered[2];
    size_t surfaceAlignment;
    int concurrentKernels;
    int ECCEnabled;
    int pciBusID;
    int pciDeviceID;
    int pciDomainID;
    int tccDriver;
    int asyncEngineCount;
    int unifiedAddressing;
    int memoryClockRate;
    int memoryBusWidth;
    int l2CacheSize;
    int maxThreadsPerMultiProcessor;
    int streamPrioritiesSupported;
    int globalL1CacheSupported;
    int localL1CacheSupported;
    size_t sharedMemPerMultiprocessor;
    int regsPerMultiprocessor;
    int managedMemSupported;
    int isMultiGpuBoard;
    int multiGpuBoardGroupID;
    int singleToDoublePrecisionPerfRatio;
    int pageableMemoryAccess;
    int concurrentManagedAccess;
    char padding[1024]; /* More room for future evolutions */
};

enum cudaComputeMode {
    cudaComputeModeDefault          = 0, 
    cudaComputeModeExclusive        = 1, 
    cudaComputeModeProhibited       = 2, 
    cudaComputeModeExclusiveProcess = 3   
};

enum cudaMemcpyKind {
    cudaMemcpyHostToHost          =   0, 
    cudaMemcpyHostToDevice        =   1, 
    cudaMemcpyDeviceToHost        =   2, 
    cudaMemcpyDeviceToDevice      =   3, 
    cudaMemcpyDefault             =   4  
};

typedef int cudaError_t;

typedef cudaError_t (*FUNPTR_cudaGetDeviceCount)(int* device_count);
typedef cudaError_t (*FUNPTR_cudaGetDeviceProperties)(
    struct cudaDeviceProp *props, int device
);
typedef cudaError_t (*FUNPTR_cudaDeviceReset)(void);
typedef cudaError_t (*FUNPTR_cudaMalloc)(void **devPtr, size_t size);
typedef cudaError_t (*FUNPTR_cudaFree)(void* devPtr);
typedef cudaError_t (*FUNPTR_cudaMemcpy)(
    void *dst, const void *src, size_t count, enum cudaMemcpyKind kind
);

/**
 * \brief Finds and initializes a function pointer to
 *  one of the functions in CUDA.
 * \details Function pointers are stored into the 
 *  CUDAContext returned by the function CUDA().
 *  If a symbol is not found, returns NL_FALSE from the
 *  calling function.
 */
#define find_cuda_func(name)                                  \
    if(                                                       \
        (                                                     \
            CUDA()->name =                                    \
            (FUNPTR_##name)nlFindFunction(                    \
		   CUDA()->DLL_cudart,#name                   \
	    )					              \
        ) == NULL                                             \
    ) {                                                       \
        nlError("nlInitExtension_CUDA: function not found", #name); \
        return NL_FALSE;                                      \
    }

/**********************************************************/
/*      CUBLAS structures and functions                   */
/**********************************************************/

struct cublasContext;
typedef struct cublasContext *cublasHandle_t;
typedef int cublasStatus_t;

typedef enum {
    CUBLAS_SIDE_LEFT =0, 
    CUBLAS_SIDE_RIGHT=1
} cublasSideMode_t; 

typedef enum {
    CUBLAS_FILL_MODE_LOWER=0, 
    CUBLAS_FILL_MODE_UPPER=1
} cublasFillMode_t;

typedef enum {
    CUBLAS_OP_N=0,  
    CUBLAS_OP_T=1,  
    CUBLAS_OP_C=2  
} cublasOperation_t;

typedef enum {
    CUBLAS_DIAG_NON_UNIT=0, 
    CUBLAS_DIAG_UNIT=1
} cublasDiagType_t; 

typedef cublasStatus_t (*FUNPTR_cublasCreate)(cublasHandle_t* handle);
typedef cublasStatus_t (*FUNPTR_cublasDestroy)(cublasHandle_t handle);

typedef cublasStatus_t (*FUNPTR_cublasGetVersion)(
    cublasHandle_t handle, int* version
);

typedef cublasStatus_t (*FUNPTR_cublasDdot)(
    cublasHandle_t handle, int n,
    const double *x, int incx,
    const double *y, int incy,
    double *result
);

typedef cublasStatus_t (*FUNPTR_cublasDcopy)(
    cublasHandle_t handle, int n,
    const double *x, int incx,
    const double *y, int incy
);

typedef cublasStatus_t (*FUNPTR_cublasDaxpy)(
    cublasHandle_t handle, int n,
    const double* alpha, 
    const double *x, int incx,
    const double *y, int incy
);

typedef cublasStatus_t (*FUNPTR_cublasDscal)(
    cublasHandle_t handle, int n,
    const double* alpha, 
    const double *x, int incx
);

typedef cublasStatus_t (*FUNPTR_cublasDnrm2)(
    cublasHandle_t handle, int n,
    const double *x, int incx,
    double* result
);
                                
typedef cublasStatus_t (*FUNPTR_cublasDdgmm)(
    cublasHandle_t handle, cublasSideMode_t mode,
    int m, int n,
    const double* A, int lda,
    const double* x, int incx,
    double* C, int ldc
);

typedef cublasStatus_t (*FUNPTR_cublasDgemv)(
    cublasHandle_t handle, 
    cublasOperation_t trans, 
    int m,
    int n,
    const double *alpha,
    const double *A,
    int lda,
    const double *x,
    int incx,
    const double *beta, 
    double *y, 
    int incy    
);

typedef cublasStatus_t (*FUNPTR_cublasDtpsv)(
    cublasHandle_t handle, cublasFillMode_t uplo,
    cublasOperation_t trans, cublasDiagType_t diag,
    int n, const double *AP,
    double* x, int incx
);


/**
 * \brief Finds and initializes a function pointer to
 *  one of the functions in CUBLAS.
 * \details Function pointers are stored into the 
 *  CUDAContext returned by the function CUDA().
 *  If a symbol is not found, returns NL_FALSE from the
 *  calling function. Here we use the functions prefixed
 *  by "_v2".
 */
#define find_cublas_func(name)	    		              \
    if(                                                       \
        (                                                     \
            CUDA()->name =                                    \
            (FUNPTR_##name)nlFindFunction(                    \
		   CUDA()->DLL_cublas,#name "_v2"             \
	    )					              \
        ) == NULL                                             \
    ) {                                                       \
        nlError("nlInitExtension_CUDA: function not found", #name); \
        return NL_FALSE;                                      \
    }

/**
 * \brief Finds and initializes a function pointer to
 *  one of the functions in CUBLAS.
 * \details Function pointers are stored into the 
 *  CUDAContext returned by the function CUDA().
 *  If a symbol is not found, returns NL_FALSE from the
 *  calling function. Here we use the functions prefixed
 *  by "_v2".
 */
#define find_cublas_func_v1(name)                             \
    if(                                                       \
        (                                                     \
            CUDA()->name =                                    \
            (FUNPTR_##name)nlFindFunction(                    \
		   CUDA()->DLL_cublas,#name                   \
	    )					              \
        ) == NULL                                             \
    ) {                                                       \
        nlError("nlInitExtension_CUDA: function not found", #name); \
        return NL_FALSE;                                      \
    }



/**********************************************************/
/*      CUSPARSE structures and functions                 */
/**********************************************************/

struct cusparseContext;
typedef struct cusparseContext *cusparseHandle_t;
typedef int cusparseStatus_t;
struct cusparseMatDescr;
typedef struct cusparseMatDescr *cusparseMatDescr_t;

typedef enum {
    CUSPARSE_MATRIX_TYPE_GENERAL = 0, 
    CUSPARSE_MATRIX_TYPE_SYMMETRIC = 1,     
    CUSPARSE_MATRIX_TYPE_HERMITIAN = 2, 
    CUSPARSE_MATRIX_TYPE_TRIANGULAR = 3 
} cusparseMatrixType_t;

typedef enum {
    CUSPARSE_INDEX_BASE_ZERO = 0, 
    CUSPARSE_INDEX_BASE_ONE = 1
} cusparseIndexBase_t;

typedef enum {
    CUSPARSE_OPERATION_NON_TRANSPOSE = 0,  
    CUSPARSE_OPERATION_TRANSPOSE = 1,  
    CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE = 2  
} cusparseOperation_t;

struct cusparseHybMat;
typedef struct cusparseHybMat *cusparseHybMat_t;

typedef enum {
    CUSPARSE_HYB_PARTITION_AUTO = 0,  
    CUSPARSE_HYB_PARTITION_USER = 1,  
    CUSPARSE_HYB_PARTITION_MAX = 2    
} cusparseHybPartition_t;

typedef cusparseStatus_t (*FUNPTR_cusparseCreate)(cusparseHandle_t* handle);
typedef cusparseStatus_t (*FUNPTR_cusparseDestroy)(cusparseHandle_t handle);
typedef cusparseStatus_t (*FUNPTR_cusparseGetVersion)(
    cusparseHandle_t handle, int* version
);
typedef cusparseStatus_t (*FUNPTR_cusparseCreateMatDescr)(
    cusparseMatDescr_t* descr
);
typedef cusparseStatus_t (*FUNPTR_cusparseDestroyMatDescr)(
    cusparseMatDescr_t descr
);
typedef cusparseStatus_t (*FUNPTR_cusparseSetMatType)(
    cusparseMatDescr_t descr, cusparseMatrixType_t mtype
);
typedef cusparseStatus_t (*FUNPTR_cusparseSetMatIndexBase)(
    cusparseMatDescr_t descr, cusparseIndexBase_t ibase
);
typedef cusparseStatus_t (*FUNPTR_cusparseDcsrmv)(
     cusparseHandle_t handle, cusparseOperation_t transA, 
     int m, int n, int nnz,
     const double *alpha, const cusparseMatDescr_t descrA, 
     const double *csrSortedValA, const int *csrSortedRowPtrA, 
     const int *csrSortedColIndA, const double *x, 
     const double *beta, double *y
);
typedef cusparseStatus_t (*FUNPTR_cusparseCreateHybMat)(
    cusparseHybMat_t *hybA
);
typedef cusparseStatus_t (*FUNPTR_cusparseDestroyHybMat)(
    cusparseHybMat_t hybA
);
typedef cusparseStatus_t (*FUNPTR_cusparseDcsr2hyb)(
    cusparseHandle_t handle,
    int m,
    int n,
    const cusparseMatDescr_t descrA,
    const double *csrSortedValA,
    const int *csrSortedRowPtrA,
    const int *csrSortedColIndA,
    cusparseHybMat_t hybA,
    int userEllWidth,
    cusparseHybPartition_t partitionType
);
typedef cusparseStatus_t (*FUNPTR_cusparseDhybmv)(
    cusparseHandle_t handle,
    cusparseOperation_t transA,
    const double *alpha,
    const cusparseMatDescr_t descrA,
    const cusparseHybMat_t hybA,
    const double *x,
    const double *beta,
    double *y
);

/**
 * \brief Finds and initializes a function pointer to
 *  one of the functions in CUSPARSE.
 * \details Function pointers are stored into the 
 *  CUDAContext returned by the function CUDA().
 *  If a symbol is not found, returns NL_FALSE from the
 *  calling function. 
 */
#define find_cusparse_func(name)                              \
    if(                                                       \
        (                                                     \
            CUDA()->name =                                    \
            (FUNPTR_##name)nlFindFunction(                    \
		   CUDA()->DLL_cusparse,#name                 \
	    )					              \
        ) == NULL                                             \
    ) {                                                       \
        nlError("nlInitExtension_CUDA : function not found", #name);	\
        return NL_FALSE;                                      \
    }


/**********************************************************/

/**
 * \brief The structure that stores the handle to 
 *  the CUDA shared object, the function pointers
 *  and the detected version.
 */
typedef struct {
    NLdll DLL_cudart;
    FUNPTR_cudaGetDeviceCount cudaGetDeviceCount;
    FUNPTR_cudaGetDeviceProperties cudaGetDeviceProperties;
    FUNPTR_cudaDeviceReset cudaDeviceReset;
    FUNPTR_cudaMalloc cudaMalloc;
    FUNPTR_cudaFree cudaFree;
    FUNPTR_cudaMemcpy cudaMemcpy;

    NLdll DLL_cublas;
    cublasHandle_t HNDL_cublas;
    FUNPTR_cublasCreate cublasCreate;
    FUNPTR_cublasDestroy cublasDestroy;
    FUNPTR_cublasGetVersion cublasGetVersion;
    FUNPTR_cublasDdot cublasDdot;
    FUNPTR_cublasDcopy cublasDcopy;
    FUNPTR_cublasDaxpy cublasDaxpy;
    FUNPTR_cublasDscal cublasDscal;    
    FUNPTR_cublasDnrm2 cublasDnrm2;
    FUNPTR_cublasDdgmm cublasDdgmm;
    FUNPTR_cublasDgemv cublasDgemv;
    FUNPTR_cublasDtpsv cublasDtpsv;
    
    NLdll DLL_cusparse;
    cusparseHandle_t HNDL_cusparse;
    FUNPTR_cusparseCreate cusparseCreate;
    FUNPTR_cusparseDestroy cusparseDestroy;
    FUNPTR_cusparseGetVersion cusparseGetVersion;
    FUNPTR_cusparseCreateMatDescr cusparseCreateMatDescr;
    FUNPTR_cusparseDestroyMatDescr cusparseDestroyMatDescr;    
    FUNPTR_cusparseSetMatType cusparseSetMatType;
    FUNPTR_cusparseSetMatIndexBase cusparseSetMatIndexBase;    
    FUNPTR_cusparseDcsrmv cusparseDcsrmv;
    FUNPTR_cusparseCreateHybMat cusparseCreateHybMat;
    FUNPTR_cusparseDestroyHybMat cusparseDestroyHybMat;
    FUNPTR_cusparseDcsr2hyb cusparseDcsr2hyb;
    FUNPTR_cusparseDhybmv cusparseDhybmv;
    
    int devID;
} CUDAContext;

/**
 * \brief Gets the CUDA context.
 * \return a pointer to the CUDA context
 */
static CUDAContext* CUDA() {
    static CUDAContext context;
    static NLboolean init = NL_FALSE;
    if(!init) {
        init = NL_TRUE;
        memset(&context, 0, sizeof(context));
    }
    return &context;
}

NLboolean nlExtensionIsInitialized_CUDA() {
    if(
	CUDA()->DLL_cudart == NULL ||
	CUDA()->cudaGetDeviceCount == NULL ||
	CUDA()->cudaGetDeviceProperties == NULL ||
	CUDA()->cudaDeviceReset == NULL ||
	CUDA()->cudaMalloc == NULL ||
	CUDA()->cudaFree == NULL ||
	CUDA()->cudaMemcpy == NULL ||
	
	CUDA()->DLL_cublas == NULL ||
	CUDA()->HNDL_cublas == NULL ||
	CUDA()->cublasCreate == NULL ||
	CUDA()->cublasDestroy == NULL ||
	CUDA()->cublasGetVersion == NULL ||
	CUDA()->cublasDdot == NULL ||
	CUDA()->cublasDcopy == NULL ||
	CUDA()->cublasDaxpy == NULL ||
	CUDA()->cublasDscal == NULL ||
	CUDA()->cublasDnrm2 == NULL ||
	CUDA()->cublasDdgmm == NULL ||
	
	CUDA()->DLL_cusparse == NULL ||
	CUDA()->HNDL_cusparse == NULL ||
	CUDA()->cusparseCreate == NULL ||
	CUDA()->cusparseDestroy == NULL ||
	CUDA()->cusparseGetVersion == NULL ||
	CUDA()->cusparseCreateMatDescr == NULL ||
	CUDA()->cusparseDestroyMatDescr == NULL ||	
	CUDA()->cusparseSetMatType == NULL ||
	CUDA()->cusparseSetMatIndexBase == NULL ||
	CUDA()->cusparseDcsrmv == NULL ||
	CUDA()->cusparseCreateHybMat == NULL ||
	CUDA()->cusparseDestroyHybMat == NULL ||
	CUDA()->cusparseDcsr2hyb == NULL ||
	CUDA()->cusparseDhybmv == NULL
    ) {
        return NL_FALSE;
    }
    return NL_TRUE;
}

static void nlTerminateExtension_CUDA(void) {
    if(!nlExtensionIsInitialized_CUDA()) {
	return;
    }

    CUDA()->cusparseDestroy(CUDA()->HNDL_cusparse);    
    nlCloseDLL(CUDA()->DLL_cusparse);
    
    CUDA()->cublasDestroy(CUDA()->HNDL_cublas);
    nlCloseDLL(CUDA()->DLL_cublas);

    CUDA()->cudaDeviceReset();    
    nlCloseDLL(CUDA()->DLL_cudart);

    memset(CUDA(), 0, sizeof(CUDAContext));
}

/**************************************************************************/

/**
 * \brief Finds the number of cores from the major and minor versions of the
 *  shader model.
 * \details Highly inspired by the helpers library in CUDA examples.
 */
static int ConvertSMVer2Cores(int major, int minor) {
    /* Defines for GPU Architecture types (using the SM version 
       to determine the # of cores per SM */
    typedef struct {
        int SM; /* 0xMm (hexadecimal notation), 
                    M = SM Major version, 
                    and m = SM minor version */
        int Cores;
    } sSMtoCores;

    sSMtoCores nGpuArchCoresPerSM[] = {
        { 0x10,  8 }, /* Tesla Generation   (SM 1.0) G80 class    */
        { 0x11,  8 }, /* Tesla Generation   (SM 1.1) G8x class    */
        { 0x12,  8 }, /* Tesla Generation   (SM 1.2) G9x class    */
        { 0x13,  8 }, /* Tesla Generation   (SM 1.3) GT200 class  */
        { 0x20, 32 }, /* Fermi Generation   (SM 2.0) GF100 class  */
        { 0x21, 48 }, /* Fermi Generation   (SM 2.1) GF10x class  */
        { 0x30, 192}, /* Kepler Generation  (SM 3.0) GK10x class  */
        { 0x35, 192}, /* Kepler Generation  (SM 3.5) GK11x class  */
	{ 0x50, 128}, /* Maxwell Generation (SM 5.0) GM10x class 
                             (yes, #cores smaller than with 3.x)  */
	{ 0x52, 128}, /* Maxwell Generation (SM 5.2) GM20x class  */
	{ 0x60, 64 }, /* Pascal Generation  (SM 6.0) GP100,GP102  
              (yes, 64, but GP100 has superfast double precision) */
	{ 0x61, 128}, /* Pascal Generation  (SM 6.1) GP104 class  
                               (but FP64 runs as 1/32 FP32 speed) */ 	
        {   -1, -1 }
    };
    int index = 0;
    while (nGpuArchCoresPerSM[index].SM != -1) {
        if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor)) {
            return nGpuArchCoresPerSM[index].Cores;
        }
        index++;
    }
    /* If we don't find the values, we default use the 
       previous one to run properly */
    nl_printf(
      "MapSMtoCores for SM %d.%d is undefined.  Default to use %d Cores/SM\n",
      major, minor, nGpuArchCoresPerSM[8].Cores
    );
    return nGpuArchCoresPerSM[8].Cores;
}

/**
 * \brief Finds among all detected GPUs the fastest one.
 * \details Highly inspired by the helpers library in
 *  CUDA examples.
 */
static int getBestDeviceID() {
    int current_device     = 0, sm_per_multiproc  = 0;
    int max_compute_perf   = 0, max_perf_device   = 0;
    int device_count       = 0, best_SM_arch      = 0;
    int compute_perf       = 0;
    struct cudaDeviceProp deviceProp;
    CUDA()->cudaGetDeviceCount(&device_count);
    /* Find the best major SM Architecture GPU device */
    while (current_device < device_count) {
        CUDA()->cudaGetDeviceProperties(&deviceProp, current_device);
        /* If this GPU is not running on Compute Mode prohibited, 
           then we can add it to the list */
        if (deviceProp.computeMode != cudaComputeModeProhibited) {
            if (deviceProp.major > 0 && deviceProp.major < 9999) {
                best_SM_arch = MAX(best_SM_arch, deviceProp.major);
            }
        }
        current_device++;
    }
    /* Find the best CUDA capable GPU device */
    current_device = 0;
    while (current_device < device_count) {
        CUDA()->cudaGetDeviceProperties(&deviceProp, current_device);
        /* If this GPU is not running on Compute Mode prohibited, 
           then we can add it to the list */
        if (deviceProp.computeMode != cudaComputeModeProhibited) {
            if (deviceProp.major == 9999 && deviceProp.minor == 9999) {
                sm_per_multiproc = 1;
            } else {
                sm_per_multiproc = ConvertSMVer2Cores(
		    deviceProp.major, deviceProp.minor
		);
            }
            compute_perf  =
		deviceProp.multiProcessorCount *
		sm_per_multiproc * deviceProp.clockRate;
            if (compute_perf  > max_compute_perf) {
                /* If we find GPU with SM major > 2, search only these */
                if (best_SM_arch > 2) {
                    /* If our device==dest_SM_arch, choose this, or else pass */
                    if (deviceProp.major == best_SM_arch) {
                        max_compute_perf  = compute_perf;
                        max_perf_device   = current_device;
                    }
                } else {
                    max_compute_perf  = compute_perf;
                    max_perf_device   = current_device;
                }
            }
        }
        ++current_device;
    }

    /**
     * Wrong ! It says 1TFlops for GTX 1080, whereas the specs say 8TFlops
     nl_printf(
	"OpenNL CUDA: maximum device single-precision Gflops=%f\n",
	(double)(2*max_compute_perf)/(double)(1e6)
    );
    */
 
    return max_perf_device;
}


/**************************************************************************/

#ifdef NL_OS_UNIX
#  define LIBPREFIX "lib"
#  ifdef NL_OS_APPLE
#      define LIBEXTENSION ".dylib"
#  else
#      define LIBEXTENSION ".so"
#  endif
#else
#  define LIBPREFIX 
#  define LIBEXTENSION ".dll"
#endif


NLboolean nlInitExtension_CUDA(void) {
    struct cudaDeviceProp deviceProp;
    int cublas_version;
    int cusparse_version;
    NLenum flags = NL_LINK_LAZY | NL_LINK_GLOBAL;
    if(nlCurrentContext == NULL || !nlCurrentContext->verbose) {
	flags |= NL_LINK_QUIET;
    }
    
    if(nlExtensionIsInitialized_CUDA()) {
	return NL_TRUE;
    }

    CUDA()->DLL_cudart = nlOpenDLL(
	LIBPREFIX "cudart" LIBEXTENSION, flags
    );

    find_cuda_func(cudaGetDeviceCount);
    find_cuda_func(cudaGetDeviceProperties);
    find_cuda_func(cudaDeviceReset);        
    find_cuda_func(cudaMalloc);
    find_cuda_func(cudaFree);
    find_cuda_func(cudaMemcpy);
    
    CUDA()->devID = getBestDeviceID();

    if(CUDA()->cudaGetDeviceProperties(&deviceProp, CUDA()->devID)) {
	nl_fprintf(stderr,"OpenNL CUDA: could not find a CUDA device\n");
	return NL_FALSE;
    }
    
    nl_printf("OpenNL CUDA: Device ID = %d\n", CUDA()->devID);
    nl_printf("OpenNL CUDA: Device name=%s\n", deviceProp.name);
    nl_printf(
	"OpenNL CUDA: Device has %d Multi-Processors, "
	"%d cores per Multi-Processor, SM %d.%d compute capabilities\n",
	deviceProp.multiProcessorCount,
	ConvertSMVer2Cores(deviceProp.major, deviceProp.minor),
	deviceProp.major, deviceProp.minor
    );

    nl_printf(
	"OpenNL CUDA: %d kB shared mem. per block, %d per MP\n",
	(int)(deviceProp.sharedMemPerBlock / 1024),
	(int)(deviceProp.sharedMemPerMultiprocessor / 1024)
    );
    
    nl_printf(
	"OpenNL CUDA: %d regs. per block, %d per MP\n",
	deviceProp.regsPerBlock,
	deviceProp.regsPerMultiprocessor	
    );

    nl_printf(
	"OpenNL CUDA: warpsize=%d\n",
	deviceProp.warpSize
    );
    
    if ((deviceProp.major * 0x10 + deviceProp.minor) < 0x11) {
        nl_fprintf(stderr, "OpenNL CUDA requires a minimum CUDA compute 1.1 capability\n");
        CUDA()->cudaDeviceReset();
	return NL_FALSE;
    }
    
    CUDA()->DLL_cublas = nlOpenDLL(
	LIBPREFIX "cublas" LIBEXTENSION, flags
    );

    find_cublas_func(cublasCreate);
    find_cublas_func(cublasDestroy);
    find_cublas_func(cublasGetVersion);    
    find_cublas_func(cublasDdot);
    find_cublas_func(cublasDaxpy);
    find_cublas_func(cublasDcopy);
    find_cublas_func(cublasDscal);
    find_cublas_func(cublasDnrm2);
    find_cublas_func(cublasDgemv);        
    find_cublas_func(cublasDtpsv);    
    find_cublas_func_v1(cublasDdgmm);

    
    if(CUDA()->cublasCreate(&CUDA()->HNDL_cublas)) {
	return NL_FALSE;
    }

    if(CUDA()->cublasGetVersion(CUDA()->HNDL_cublas, &cublas_version)) {
	return NL_FALSE;
    }
    nl_printf("OpenNL CUDA: cublas version = %d\n", cublas_version);
    
    CUDA()->DLL_cusparse = nlOpenDLL(
	LIBPREFIX "cusparse" LIBEXTENSION, flags
    );
    find_cusparse_func(cusparseCreate);
    find_cusparse_func(cusparseDestroy);
    find_cusparse_func(cusparseGetVersion);
    find_cusparse_func(cusparseCreateMatDescr);
    find_cusparse_func(cusparseDestroyMatDescr);    
    find_cusparse_func(cusparseSetMatType);
    find_cusparse_func(cusparseSetMatIndexBase);
    find_cusparse_func(cusparseDcsrmv);
    find_cusparse_func(cusparseCreateHybMat);
    find_cusparse_func(cusparseDestroyHybMat);
    find_cusparse_func(cusparseDcsr2hyb);
    find_cusparse_func(cusparseDhybmv);                    
    
    if(CUDA()->cusparseCreate(&CUDA()->HNDL_cusparse)) {
	return NL_FALSE;
    }
    if(CUDA()->cusparseGetVersion(CUDA()->HNDL_cusparse, &cusparse_version)) {
	return NL_FALSE;
    }
    nl_printf("OpenNL CUDA: cusparse version = %d\n", cusparse_version);
    
    if(!nlExtensionIsInitialized_CUDA()) {
	return NL_FALSE;
    }

    atexit(nlTerminateExtension_CUDA);
    return NL_TRUE;
    
}

static void nlCUDACheckImpl(int status, int line) {
    if(status != 0) {
	nl_fprintf(stderr,"nl_cuda.c:%d fatal error %d\n",line, status);
	CUDA()->cudaDeviceReset();    	
	exit(-1);
    }
}

#define nlCUDACheck(status) nlCUDACheckImpl(status, __LINE__)

/**************************************************************************/

/**
 * Abstract matrix interface for a CRS matrix stored on the GPU.
 */
typedef struct {
    NLuint m;
    NLuint n;
    NLenum type;
    NLDestroyMatrixFunc destroy_func;
    NLMultMatrixVectorFunc mult_func;
    cusparseMatDescr_t descr;
    NLuint nnz;
    int* colind;
    int* rowptr;
    double* val;
    cusparseHybMat_t hyb;
} NLCUDASparseMatrix;


/**
 * \brief Deallocates just the CRS part of a CUDA matrix.
 */
static void nlCRSMatrixCUDADestroyCRS(NLCUDASparseMatrix* Mcuda) {
    if(!nlExtensionIsInitialized_CUDA()) {
	return;
    }
    if(Mcuda->colind != NULL) {
	nlCUDACheck(CUDA()->cudaFree(Mcuda->colind));
	Mcuda->colind = NULL;
    }
    if(Mcuda->rowptr != NULL) {
	nlCUDACheck(CUDA()->cudaFree(Mcuda->rowptr));
	Mcuda->rowptr = NULL;
    }
    if(Mcuda->val != NULL) {
	nlCUDACheck(CUDA()->cudaFree(Mcuda->val));
	Mcuda->val = NULL;
    }
}

static void nlCRSMatrixCUDADestroy(NLCUDASparseMatrix* Mcuda) {
    if(!nlExtensionIsInitialized_CUDA()) {
	return;
    }
    if(Mcuda->hyb != NULL) {
	nlCUDACheck(CUDA()->cusparseDestroyHybMat(Mcuda->hyb));
    }
    nlCRSMatrixCUDADestroyCRS(Mcuda);
    nlCUDACheck(CUDA()->cusparseDestroyMatDescr(Mcuda->descr));
    memset(Mcuda, 0, sizeof(*Mcuda));
}

static void nlCRSMatrixCUDAMult(
    NLCUDASparseMatrix* Mcuda, const double* x, double* y
) {
    const double one = 1;
    const double zero = 0;
    if(Mcuda->hyb != NULL) {
	nlCUDACheck(
	    CUDA()->cusparseDhybmv(
		CUDA()->HNDL_cusparse,
		CUSPARSE_OPERATION_NON_TRANSPOSE,
		&one,
		Mcuda->descr,
		Mcuda->hyb,
		x,
		&zero,
		y
	    )
	);
    } else {
	nlCUDACheck(
	    CUDA()->cusparseDcsrmv(
		CUDA()->HNDL_cusparse,
		CUSPARSE_OPERATION_NON_TRANSPOSE,
		(int)Mcuda->m,
		(int)Mcuda->n,
		(int)Mcuda->nnz,
		&one,
		Mcuda->descr,
		Mcuda->val,
		Mcuda->rowptr,
		Mcuda->colind,
		x,
		&zero,
		y
		)
	    );
    }
    nlCUDABlas()->flops += (NLulong)(2*Mcuda->nnz);
}

NLMatrix nlCUDAMatrixNewFromCRSMatrix(NLMatrix M_in) {
    NLCUDASparseMatrix* Mcuda = NL_NEW(NLCUDASparseMatrix);
    NLCRSMatrix* M = (NLCRSMatrix*)(M_in);
    size_t colind_sz, rowptr_sz, val_sz;
    nl_assert(M_in->type == NL_MATRIX_CRS);
    nlCUDACheck(CUDA()->cusparseCreateMatDescr(&Mcuda->descr));
    if(M->symmetric_storage) {
	nlCUDACheck(CUDA()->cusparseSetMatType(
			Mcuda->descr, CUSPARSE_MATRIX_TYPE_SYMMETRIC)
	);
    } else {
	nlCUDACheck(CUDA()->cusparseSetMatType(
			Mcuda->descr, CUSPARSE_MATRIX_TYPE_GENERAL)
	);	
    }
    nlCUDACheck(CUDA()->cusparseSetMatIndexBase(
		    Mcuda->descr, CUSPARSE_INDEX_BASE_ZERO)
    );	
    Mcuda->m = M->m;
    Mcuda->n = M->n;
    Mcuda->nnz = nlCRSMatrixNNZ(M);

    colind_sz = (size_t)Mcuda->nnz*sizeof(int);
    rowptr_sz = (size_t)(Mcuda->m+1)*sizeof(int);
    val_sz    = (size_t)Mcuda->nnz*sizeof(double);

    nlCUDACheck(CUDA()->cudaMalloc((void**)&Mcuda->colind,colind_sz));
    nlCUDACheck(CUDA()->cudaMalloc((void**)&Mcuda->rowptr,rowptr_sz));
    nlCUDACheck(CUDA()->cudaMalloc((void**)&Mcuda->val,val_sz));
    nlCUDACheck(CUDA()->cudaMemcpy(
      Mcuda->colind, M->colind, colind_sz, cudaMemcpyHostToDevice)
    );
    nlCUDACheck(CUDA()->cudaMemcpy(
      Mcuda->rowptr, M->rowptr, rowptr_sz, cudaMemcpyHostToDevice)
    );
    nlCUDACheck(CUDA()->cudaMemcpy(
      Mcuda->val, M->val, val_sz, cudaMemcpyHostToDevice)
    );
    Mcuda->hyb=NULL;
    if(!M->symmetric_storage) {
	nlCUDACheck(CUDA()->cusparseCreateHybMat(&Mcuda->hyb));
	nlCUDACheck(CUDA()->cusparseDcsr2hyb(
			CUDA()->HNDL_cusparse,
			(int)M->m,
			(int)M->n,
			Mcuda->descr,
			Mcuda->val,
			Mcuda->rowptr,
			Mcuda->colind,
			Mcuda->hyb,
			0,
			CUSPARSE_HYB_PARTITION_AUTO 
	));
	/* We no longer need the CRS part */
	nlCRSMatrixCUDADestroyCRS(Mcuda);	
    }
    Mcuda->type=NL_MATRIX_OTHER;
    Mcuda->destroy_func=(NLDestroyMatrixFunc)nlCRSMatrixCUDADestroy;
    Mcuda->mult_func=(NLMultMatrixVectorFunc)nlCRSMatrixCUDAMult;
    return (NLMatrix)Mcuda;
}

/**************************************************************************/

/**
 * Abstract matrix interface for a diagonal matrix stored on the GPU.
 */
typedef struct {
    NLuint m;
    NLuint n;
    NLenum type;
    NLDestroyMatrixFunc destroy_func;
    NLMultMatrixVectorFunc mult_func;
    double* val;
} NLDiagonalMatrixCUDA;

static void nlDiagonalMatrixCUDADestroy(NLDiagonalMatrixCUDA* Mcuda) {
    if(!nlExtensionIsInitialized_CUDA()) {
	return;
    }
    nlCUDACheck(CUDA()->cudaFree(Mcuda->val));
    memset(Mcuda, 0, sizeof(*Mcuda));
}

static void nlDiagonalMatrixCUDAMult(
    NLDiagonalMatrixCUDA* Mcuda, const double* x, double* y
) {
    int N = (int)Mcuda->n;
    /*
     * vector x vector component-wise product implemented
     * using diagonal matrix x matrix function.
     */
    nlCUDACheck(CUDA()->cublasDdgmm(
	CUDA()->HNDL_cublas, CUBLAS_SIDE_LEFT,
	N, 1,
	x, N,
	Mcuda->val, 1,
	y, N
    ));
    nlCUDABlas()->flops += (NLulong)N;
}

static NLMatrix nlDiagonalMatrixCUDANew(const double* diag, NLuint n) {
    NLDiagonalMatrixCUDA* Mcuda = NL_NEW(NLDiagonalMatrixCUDA);
    Mcuda->m = n;
    Mcuda->n = n;
    Mcuda->type = NL_MATRIX_OTHER;
    nlCUDACheck(CUDA()->cudaMalloc(
       (void**)&Mcuda->val, n*sizeof(double))
    );
    nlCUDACheck(CUDA()->cudaMemcpy(
       Mcuda->val, diag, n*sizeof(double), cudaMemcpyHostToDevice)
    );
    Mcuda->destroy_func=(NLDestroyMatrixFunc)nlDiagonalMatrixCUDADestroy;
    Mcuda->mult_func=(NLMultMatrixVectorFunc)nlDiagonalMatrixCUDAMult;
    return (NLMatrix)Mcuda;
}

NLMatrix nlCUDAJacobiPreconditionerNewFromCRSMatrix(NLMatrix M_in) {
    NLuint N = M_in->n;
    NLuint i,jj;
    double* diag = NULL;
    NLMatrix result = NULL;
    NLCRSMatrix* M = (NLCRSMatrix*)(M_in);
    nl_assert(M_in->type == NL_MATRIX_CRS);
    diag = NL_NEW_ARRAY(double,N);
    for(i=0; i<N; ++i) {
	for(jj=M->rowptr[i]; jj<M->rowptr[i+1]; ++jj) {
	    if(M->colind[jj] == i) {
		diag[i] = M->val[jj];
	    }
	}
    }
    for(i=0; i<N; ++i) {
	diag[i] = ((diag[i] == 0.0) ? 1.0 : 1.0 / diag[i]);
    }
    result = nlDiagonalMatrixCUDANew(diag, N);
    NL_DELETE_ARRAY(diag);
    return result;
}

/**************************************************************************/

static void* cuda_blas_malloc(
    NLBlas_t blas, NLmemoryType type, size_t size
) {
    void* result = NULL;
    blas->used_ram[type] += (NLulong)size;
    blas->max_used_ram[type] = MAX(
	blas->max_used_ram[type],blas->used_ram[type]
    );
    if(type == NL_HOST_MEMORY) {
	result = malloc(size);
    } else {
	nlCUDACheck(CUDA()->cudaMalloc(&result,size));	
    }
    return result;
}

static void cuda_blas_free(
    NLBlas_t blas, NLmemoryType type, size_t size, void* ptr
) {
    blas->used_ram[type] -= (NLulong)size;
    if(type == NL_HOST_MEMORY) {
	free(ptr);
    } else {
	nlCUDACheck(CUDA()->cudaFree(ptr));	
    }
}

static void cuda_blas_memcpy(
    NLBlas_t blas,
    void* to, NLmemoryType to_type,
    void* from, NLmemoryType from_type,
    size_t size
) {
    enum cudaMemcpyKind kind = cudaMemcpyDefault;
    nl_arg_used(blas);
    if(from_type == NL_HOST_MEMORY) {
	if(to_type == NL_HOST_MEMORY) {
	    kind = cudaMemcpyHostToHost;
	} else {
	    kind = cudaMemcpyHostToDevice;	    
	}
    } else {
	if(to_type == NL_HOST_MEMORY) {
	    kind = cudaMemcpyDeviceToHost;
	} else {
	    kind = cudaMemcpyDeviceToDevice;	    
	}
    }
    nlCUDACheck(CUDA()->cudaMemcpy(to, from, size, kind));
}

static void cuda_blas_dcopy(
    NLBlas_t blas, int n, const double *x, int incx, double *y, int incy    
) {
    nl_arg_used(blas);
    CUDA()->cublasDcopy(CUDA()->HNDL_cublas,n,x,incx,y,incy);
}

static double cuda_blas_ddot(
    NLBlas_t blas, int n, const double *x, int incx, const double *y, int incy
) {
    double result = 0.0;
    blas->flops += (NLulong)(2*n);
    CUDA()->cublasDdot(CUDA()->HNDL_cublas,n,x,incx,y,incy,&result);
    return result;
}

static double cuda_blas_dnrm2(
    NLBlas_t blas, int n, const double *x, int incx
) {
    double result = 0.0;
    blas->flops += (NLulong)(2*n);
    CUDA()->cublasDnrm2(CUDA()->HNDL_cublas,n,x,incx,&result);    
    return result;
}

static void cuda_blas_daxpy(
    NLBlas_t blas, int n,
    double a, const double *x, int incx, double *y, int incy
) {
    blas->flops += (NLulong)(2*n);
    CUDA()->cublasDaxpy(CUDA()->HNDL_cublas,n,&a,x,incx,y,incy);        
}

static void cuda_blas_dscal(
    NLBlas_t blas, int n, double a, double *x, int incx    
) {
    blas->flops += (NLulong)n;
    CUDA()->cublasDscal(CUDA()->HNDL_cublas,n,&a,x,incx);            
}


static void cuda_blas_dgemv(
    NLBlas_t blas, MatrixTranspose trans, int m, int n, double alpha,
    const double *A, int ldA, const double *x, int incx,
    double beta, double *y, int incy 
) {
    nl_arg_used(blas);
    /* TODO: update FLOPS */
    CUDA()->cublasDgemv(
	CUDA()->HNDL_cublas, (cublasOperation_t)trans,
	m, n, &alpha, A, ldA, x, incx, &beta, y, incy
    );
}

static void cuda_blas_dtpsv(
    NLBlas_t blas, MatrixTriangle uplo, MatrixTranspose trans,
    MatrixUnitTriangular diag, int n, const double *AP,
    double *x, int incx 
) {
    nl_arg_used(blas);
    /* TODO: update FLOPS */
    CUDA()->cublasDtpsv(
	CUDA()->HNDL_cublas,
	(cublasFillMode_t)uplo,
	(cublasOperation_t)trans,
	(cublasDiagType_t)diag, n,
	AP, x, incx	
    );
}


NLBlas_t nlCUDABlas() {
    static NLboolean initialized = NL_FALSE;
    static struct NLBlas blas;
    if(!initialized) {
	memset(&blas, 0, sizeof(blas));
	blas.has_unified_memory = NL_FALSE;
	blas.Malloc = cuda_blas_malloc;
	blas.Free = cuda_blas_free;
	blas.Memcpy = cuda_blas_memcpy;
	blas.Dcopy = cuda_blas_dcopy;
	blas.Ddot = cuda_blas_ddot;
	blas.Dnrm2 = cuda_blas_dnrm2;
	blas.Daxpy = cuda_blas_daxpy;
	blas.Dscal = cuda_blas_dscal;
	blas.Dgemv = cuda_blas_dgemv;
	blas.Dtpsv = cuda_blas_dtpsv;
	nlBlasResetStats(&blas);
	initialized = NL_TRUE;
    }
    return &blas;
}


/**************************************************************************/
