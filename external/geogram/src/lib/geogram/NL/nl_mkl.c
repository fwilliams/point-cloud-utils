/*
 *  Copyright (c) 2004-2010, Bruno Levy
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

#include "nl_mkl.h"
#include "nl_context.h"

/**
 * \file nl_mkl.c
 * \brief Weak-coupling adapter to call MKL from OpenNL.
 */

typedef unsigned int MKL_INT;

typedef void (*FUNPTR_mkl_cspblas_dcsrgemv)(
    const char *transa, const MKL_INT *m, const double *a,
    const MKL_INT *ia, const MKL_INT *ja, const double *x, double *y
);

typedef void (*FUNPTR_mkl_cspblas_dcsrsymv)(
    const char *transa, const MKL_INT *m, const double *a,
    const MKL_INT *ia, const MKL_INT *ja, const double *x, double *y
 );

/**
 * \brief The structure that stores the handle to 
 *  the MKL shared object, the function pointers
 *  and the detected version.
 */
typedef struct {
    NLdll DLL_mkl_intel_lp64;
    NLdll DLL_mkl_intel_thread;
    NLdll DLL_mkl_core;
    NLdll DLL_iomp5;

    FUNPTR_mkl_cspblas_dcsrgemv mkl_cspblas_dcsrgemv;
    FUNPTR_mkl_cspblas_dcsrsymv mkl_cspblas_dcsrsymv;   
} MKLContext;

/**
 * \brief Gets the MKL context.
 * \return a pointer to the MKL context
 */
static MKLContext* MKL() {
    static MKLContext context;
    static NLboolean init = NL_FALSE;
    if(!init) {
        init = NL_TRUE;
        memset(&context, 0, sizeof(context));
    }
    return &context;
}

NLboolean nlExtensionIsInitialized_MKL() {
    if(
	MKL()->DLL_iomp5 == NULL ||
	MKL()->DLL_mkl_core == NULL ||
	MKL()->DLL_mkl_intel_thread == NULL ||
	MKL()->DLL_mkl_intel_lp64 == NULL ||
	MKL()->mkl_cspblas_dcsrgemv == NULL ||
	MKL()->mkl_cspblas_dcsrsymv == NULL 	
    ) {
        return NL_FALSE;
    }
    return NL_TRUE;
}

/**
 * \brief Finds and initializes a function pointer to
 *  one of the functions in MKL.
 * \details Function pointers are stored into the 
 *  MKLContext returned by the function MKL().
 *  If a symbol is not found, returns NL_FALSE from the
 *  calling function.
 */
#define find_mkl_func(name)                                  \
    if(                                                      \
        (                                                    \
            MKL()->name =                                    \
            (FUNPTR_##name)nlFindFunction(                   \
		   MKL()->DLL_mkl_intel_lp64,#name           \
	    )					             \
        ) == NULL                                            \
    ) {                                                      \
        nlError("nlInitExtension_MKL","function not found"); \
        return NL_FALSE;                                     \
    }

static void nlTerminateExtension_MKL(void) {
    if(!nlExtensionIsInitialized_MKL()) {
	return;
    }
    nlCloseDLL(MKL()->DLL_mkl_intel_lp64);
    nlCloseDLL(MKL()->DLL_mkl_intel_thread);
    nlCloseDLL(MKL()->DLL_mkl_core);
    nlCloseDLL(MKL()->DLL_iomp5);
    memset(MKL(), 0, sizeof(MKLContext));
    
}

NLMultMatrixVectorFunc NLMultMatrixVector_MKL = NULL;

static void NLMultMatrixVector_MKL_impl(NLMatrix M_in, const double* x, double* y) {
    NLCRSMatrix* M = (NLCRSMatrix*)(M_in);
    nl_debug_assert(M_in->type == NL_MATRIX_CRS);
    if(M->symmetric_storage) {
	MKL()->mkl_cspblas_dcsrsymv(
	    "N", /* No transpose */
	    &M->m,
	    M->val,
	    M->rowptr,
	    M->colind,
	    x,
	    y
	);
    } else {
	MKL()->mkl_cspblas_dcsrgemv(
	    "N", /* No transpose */
	    &M->m,
	    M->val,
	    M->rowptr,
	    M->colind,
	    x,
	    y
	);
    }
}


#define INTEL_PREFIX "/opt/intel/"
#define LIB_DIR "lib/intel64/"
#define MKL_PREFIX  INTEL_PREFIX "mkl/" LIB_DIR

NLboolean nlInitExtension_MKL(void) {
    NLenum flags = NL_LINK_LAZY | NL_LINK_GLOBAL;
    if(nlCurrentContext == NULL || !nlCurrentContext->verbose) {
	flags |= NL_LINK_QUIET;
    }
    
    if(MKL()->DLL_mkl_intel_lp64 != NULL) {
        return nlExtensionIsInitialized_MKL();
    }
    
    MKL()->DLL_iomp5 = nlOpenDLL(
	INTEL_PREFIX LIB_DIR "libiomp5.so",
	flags
    );    
    MKL()->DLL_mkl_core = nlOpenDLL(
	MKL_PREFIX "libmkl_core.so",
	flags
    );    
    MKL()->DLL_mkl_intel_thread = nlOpenDLL(
	MKL_PREFIX "libmkl_intel_thread.so",
	flags
    );    
    MKL()->DLL_mkl_intel_lp64 = nlOpenDLL(
	MKL_PREFIX "libmkl_intel_lp64.so",
	flags
    );
    
    if(
	MKL()->DLL_iomp5 == NULL ||
	MKL()->DLL_mkl_core == NULL ||
	MKL()->DLL_mkl_intel_thread == NULL ||
	MKL()->DLL_mkl_intel_lp64 == NULL
    ) {
        return NL_FALSE;
    }

    find_mkl_func(mkl_cspblas_dcsrgemv);
    find_mkl_func(mkl_cspblas_dcsrsymv);

    if(nlExtensionIsInitialized_MKL()) {
	NLMultMatrixVector_MKL = NLMultMatrixVector_MKL_impl;
    }
    
    atexit(nlTerminateExtension_MKL);
    return NL_TRUE;
}

/*************************************************************************/

