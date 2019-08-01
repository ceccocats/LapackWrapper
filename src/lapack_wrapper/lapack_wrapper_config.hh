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

///
/// file: lapack_wrapper_config.hh
///

#ifndef LAPACK_WRAPPER_CONFIG_HH
#define LAPACK_WRAPPER_CONFIG_HH

/*\
 ! values
 !  LAPACK_WRAPPER_USE_ACCELERATE, LAPACK_WRAPPER_USE_ATLAS, LAPACK_WRAPPER_USE_OPENBLAS, 
 !  LAPACK_WRAPPER_USE_LAPACK, LAPACK_WRAPPER_USE_MKL
\*/
#define LAPACK_WRAPPER_USE_ACCELERATE 1

/*\
 ! values
 !  LAPACK_WRAPPER_USE_THREAD, LAPACK_WRAPPER_DO_NOT_USE_CXX11
\*/
#define LAPACK_WRAPPER_USE_THREAD 1

/*\
 ! define
 ! LAPACK_WRAPPER_DO_NOT_USE_SYSTEM_OPENBLAS
 ! if you want to use locally compiled openblas
\*/
//#define LAPACK_WRAPPER_DO_NOT_USE_SYSTEM_OPENBLAS 1

// select computer architecture
#if defined(__APPLE__) && defined(__MACH__)
  // osx architecture
  #define LAPACK_WRAPPER_OS_OSX 1
  #if defined(__i386__)
    #define LAPACK_WRAPPER_ARCH32 1
  #elif defined(__x86_64__)
    #define LAPACK_WRAPPER_ARCH64 1
  #endif
#elif defined(__unix__)
  // linux architecture
  #define LAPACK_WRAPPER_OS_LINUX 1
  #if defined(__i386__)
    #define LAPACK_WRAPPER_ARCH32 1
  #elif defined(__x86_64__)
    #define LAPACK_WRAPPER_ARCH64 1
  #endif
#elif defined(_WIN32) || defined(WIN32) || defined(_WIN64) || defined(WIN64)
  // windows architecture
  #define LAPACK_WRAPPER_OS_WINDOWS 1
  #if defined(_M_X64) || defined(_M_AMD64)
    #define LAPACK_WRAPPER_ARCH64 1
  #else
    #define LAPACK_WRAPPER_ARCH32 1
  #endif
  #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
  #endif
  #include <windows.h>
#else
  #error "unsupported OS!"
#endif

// check if compiler is C++11
#if (defined(_MSC_VER) &&  _MSC_VER >= 1800) || \
    (defined(__cplusplus) && __cplusplus > 199711L)
  #ifndef LAPACK_WRAPPER_DO_NOT_USE_CXX11
    #define LAPACK_WRAPPER_USE_CXX11
  #endif
#else
  // "LapackWrapper libray compiled without c++11 support, cannot use thread"
  // not C++11 compiler
  #ifndef nullptr
    #define nullptr NULL
  #endif
#endif

#define LAPACK_WRAPPER_PURE_VIRTUAL = 0
#if defined(LAPACK_WRAPPER_USE_CXX11) && !defined(LAPACK_WRAPPER_OS_WINDOWS)
  #define LAPACK_WRAPPER_OVERRIDE  override
  #define LAPACK_WRAPPER_CONSTEXPR constexpr
  #ifdef __clang__
    #pragma clang diagnostic ignored "-Wc++98-compat"
  #endif
#else
  #define LAPACK_WRAPPER_OVERRIDE
  #define LAPACK_WRAPPER_CONSTEXPR
#endif

#endif

///
/// eof: lapack_wrapper_config.hh
///
