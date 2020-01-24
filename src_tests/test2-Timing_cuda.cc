/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  Copyright (C) 2017                                                      |
 |                                                                          |
 |         , __                 , __                                        |
 |        /|/  \               /|/  \                                       |
 |         | __/ _   ,_         | __/ _   ,_                                |
 |         |   \|/  /  |  |   | |   \|/  /  |  |   |                        |
 |         |(__/|__/   |_/ \_/|/|(__/|__/   |_/ \_/|/                       |
 |                           /|                   /|                        |
 |                           \|                   \|                        |
 |                                                                          |
 |      Enrico Bertolazzi                                                   |
 |      Dipartimento di Ingegneria Industriale                              |
 |      Universita` degli Studi di Trento                                   |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/


#include <iostream>
#include <vector>
#include <random>
#include <lapack_wrapper/lapack_wrapper.hh>
#include <lapack_wrapper/lapack_wrapper++.hh>
#include <lapack_wrapper/TicToc.hh>

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wmissing-noreturn"
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic ignored "-Wdeprecated"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wc++98-compat"
#pragma clang diagnostic ignored "-Wc++98-compat-pedantic"
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wglobal-constructors"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#pragma clang diagnostic ignored "-Wold-style-cast"
#pragma clang diagnostic ignored "-Wmissing-noreturn"
#pragma clang diagnostic ignored "-Wfloat-equal"
#pragma clang diagnostic ignored "-Wdocumentation"
#pragma clang diagnostic ignored "-Wdocumentation-deprecated-sync"
#pragma clang diagnostic ignored "-Wused-but-marked-unused"
#pragma clang diagnostic ignored "-Wdeprecated"
#pragma clang diagnostic ignored "-Wundefined-func-template"
#pragma clang diagnostic ignored "-Wextra-semi"
#endif

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wzero-as-null-pointer-constant"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wcast-qual"
#endif
#ifdef __clang__
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wc99-extensions"
#pragma clang diagnostic ignored "-Wzero-as-null-pointer-constant"
#pragma clang diagnostic ignored "-Wreserved-id-macro"
#pragma clang diagnostic ignored "-Wshadow"
#pragma clang diagnostic ignored "-Wunused-template"
#pragma clang diagnostic ignored "-Wcast-qual"
#endif

// cuda includes
#include <cuda_runtime.h>
#include "cublas_v2.h"


using namespace std;
typedef double real_type;

using lapack_wrapper::integer;

static unsigned seed1 = 2;
static std::mt19937 generator(seed1);

static
real_type
rand( real_type xmin, real_type xmax ) {
  real_type random = real_type(generator())/generator.max();
  return xmin + (xmax-xmin)*random;
}

using namespace lapack_wrapper;

#define N_TIMES 1000000

template <int N>
void
testN() {

  cout << "\nSize N = " << N << "\n" << flush;

  real_type *M1, *M2, *M3;          // host pointers
  real_type *M1_d, *M2_d, *M3_d;    // device pointers

  // alloc once
  int sizeBytes = N*N*sizeof(real_type);
  cudaHostAlloc((void **)&M1,  sizeBytes,  cudaHostAllocMapped);
  cudaHostAlloc((void **)&M2,  sizeBytes,  cudaHostAllocMapped);
  cudaHostAlloc((void **)&M3,  sizeBytes,  cudaHostAllocMapped);

  // device pointer from host memory. No allocation or memcpy
  cudaHostGetDevicePointer((void **)&M1_d,  (void *) M1, 0);
  cudaHostGetDevicePointer((void **)&M2_d,  (void *) M2, 0);
  cudaHostGetDevicePointer((void **)&M3_d,  (void *) M3, 0);

  // fill random on CPU
  for ( int i = 0; i < N; ++i ) {
    for ( int j = 0; j < N; ++j ) {
      M1[i+j*N] = rand(-1,1);
      M2[i+j*N] = rand(-1,1);
      M3[i+j*N] = rand(-1,1);
    }
  }

  TicToc tm;


  // ===========================================================================
  // GEMM on cpu
  tm.tic();
  for ( int i = 0; i < N_TIMES; ++i ) {
    gemm(
      NO_TRANSPOSE, NO_TRANSPOSE,
      N, N, N,
      -1.0, M1, N,
      M2, N,
      1.0, M3, N
    );
    copy( N*N, M3, 1, M2, 1);
  }
  tm.toc();
  cout << "CPU MULT = " << tm.elapsed_ms() << " [ms] (lapack)\n";

  // ===========================================================================
  // GEMM on gpu

  // init cublas
  cublasStatus_t stat;
  cublasHandle_t handle;
  stat = cublasCreate(&handle);
  if (stat != CUBLAS_STATUS_SUCCESS) {
    printf ("CUBLAS initialization failed\n");
  }

  tm.tic();
  for ( int i = 0; i < N_TIMES; ++i ) {
    double alpha = -1.0;
    double beta = 1.0;
    cublasDgemm(handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      N, N, N,
      &alpha, M1_d, N,
      M2_d, N,
      &beta, M3_d, N
    );
    copy( N*N, M3, 1, M2, 1);
  }
  tm.toc();
  cout << "GPU MULT = " << tm.elapsed_ms() << " [ms] (lapack)\n";

  // deallocate
  cublasDestroy(handle);
  cudaFree(M1);
  cudaFree(M2);
  cudaFree(M3);

  cout << "All done!\n" << flush;
}



int
main() {

  testN<2>();
  testN<3>();
  testN<4>();
  testN<5>();
  testN<6>();
  testN<7>();
  testN<8>();
  testN<16>();
  testN<32>();
  testN<64>();

  cout << "\n\nAll done!\n" << flush;

  return 0;
}
