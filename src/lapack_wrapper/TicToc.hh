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
 |      Via Sommarive 9, I-38123 Povo, Trento, Italy                        |
 |      email: enrico.bertolazzi@unitn.it                                   |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifndef TIC_TOC_HH
#define TIC_TOC_HH

#include "../lapack_wrapper_config.hh"

#ifdef LAPACK_WRAPPER_OS_WINDOWS
  #include <windows.h>
  class TicToc {

    typedef double real_type;
    LARGE_INTEGER frequency; // ticks per second
    LARGE_INTEGER t1, t2;    // ticks
    real_type elapsed_time;

    TicToc( TicToc const & );
    TicToc const & operator = ( TicToc const & ) const;

  public:

    TicToc()
    : elapsed_time(0)
    { QueryPerformanceFrequency(&frequency); tic(); }

    ~TicToc() {}

    void
    tic()
    { QueryPerformanceCounter(&t1); }

    void
    toc() {
      QueryPerformanceCounter(&t2);
      elapsed_time = (t2.QuadPart - t1.QuadPart) * 1000.0 / frequency.QuadPart;;
    }

    real_type
    elapsed_s() const
    { return 1e-3*elapsed_time; }

    real_type
    elapsed_ms() const
    { return elapsed_time; }

  };

  inline
  void
  sleep_for_seconds( unsigned s )
  { Sleep(DWORD(s)*1000); }

  inline
  void
  sleep_for_milliseconds( unsigned ms )
  { Sleep(DWORD(ms)); }

#else

  #include <chrono>
  #include <thread>

  class TicToc {

    typedef double real_type;
    #ifdef TIC_TOC_USE_HRC
    typedef std::chrono::high_resolution_clock clock;
    #else
    typedef std::chrono::steady_clock clock;
    #endif

    using elapsed_resolution = std::chrono::microseconds;

    clock::time_point start_time;
    clock::time_point stop_time;

    elapsed_resolution elapsed_time;

    TicToc( TicToc const & );
    TicToc const & operator = ( TicToc const & ) const;

   public:

    TicToc()
    : elapsed_time(0)
    { tic(); }

    ~TicToc() {}

    void
    tic()
    { start_time = clock::now(); }

    void
    toc() {
      stop_time    = clock::now();
      elapsed_time = std::chrono::duration_cast<elapsed_resolution>(stop_time - start_time);
    }

    real_type
    elapsed_s() const
    { return 1e-6*elapsed_time.count(); }

    real_type
    elapsed_ms() const
    { return 1e-3*elapsed_time.count(); }
  };

  inline
  void
  sleep_for_seconds( unsigned s )
  { std::this_thread::sleep_for(std::chrono::seconds(s)); }

  inline
  void
  sleep_for_milliseconds( unsigned ms )
  { std::this_thread::sleep_for(std::chrono::milliseconds(ms)); }

#endif

#endif
