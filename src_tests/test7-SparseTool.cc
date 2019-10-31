/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  SparseTool   : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH UMFPACK  |
 |                                                                          |
 |  date         : 2008, 14 April                                           |
 |  version      : 1.0                                                      |
 |  file         : test.UMF.cc                                              |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Dipartimento di Ingegneria Industriale                   |
 |                 Universita` degli Studi di Trento                        |
 |                 email : enrico.bertolazzi@unitn.it                       |
 |                                                                          |
 |  purpose:                                                                |
 |                                                                          |
 |    Test the interface with UMFPACK reading a matrix from a MatrixMarket  |
 |    file and solving a linear system.                                     |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#define SPARSETOOL_DEBUG
#include <sparse_tool/sparse_tool.hh>
#include <sparse_tool/sparse_tool_extra.hh>
#include <sparse_tool/sparse_tool_iterative.hh>
#include <sparse_tool/sparse_tool_matrix_market.hh>

#include <lapack_wrapper/TicToc.hh>

#include <sparse_tool/interfaces/MA41.hh>
#include <sparse_tool/interfaces/MA48.hh>
#include <sparse_tool/interfaces/SuperLU.hh>
#include <sparse_tool/interfaces/mkl_pardiso.hh>
#include <sparse_tool/interfaces/UMF.hh>

#include <fstream>
#include <iostream>

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshorten-64-to-32"
#endif

#include "izstream.hh"

#ifdef __clang__
#pragma clang diagnostic pop
#endif

using namespace SparseToolLoad;
using namespace std;
using namespace zstream;

void
testSparseTool( istream & mm_file ) {
  TicToc                     tm;
  MatrixMarket               mm;
  ILDUpreconditioner<double> ildu;
  //ILDUKpreconditioner<double> ildu;
  CCoorMatrix<double>        A;
  CRowMatrix<double>         B;
  SparsePattern              SP;
  Vector<double>             x, rhs, exact, resid;

  fmt::print("read matrix...");
  mm.read( mm_file, A );
  fmt::print("done\n{}",mm);
  
  SP . resize( A );
  B  . resize( SP );
  //Spy( mm_file + ".eps" , A, 15.0 );
  
  exact . resize( A.numRows() );
  x     . resize( A.numRows() );
  rhs   . resize( A.numRows() );
  resid . resize( A.numRows() );

  exact = 1;
  rhs   = A * exact;

#if 1

  x = 0;
  fmt::print("factorize (ildu) ...");
  tm.tic();
  ildu.build(A);
  tm.toc();
  fmt::print(" {} [ms] done\n", tm.elapsed_ms());
  fmt::print("solve (bicgstab) ... ");
  tm.tic();

  double   epsi       = 1e-15;
  unsigned maxIter    = 200;
  //unsigned maxSubIter = 50;
  unsigned iter;
  double   res = bicgstab( A, rhs, x, ildu, epsi, maxIter, iter, &cout );
  //double   res = gmres( A, rhs, x, ildu, epsi, maxSubIter, maxIter, iter, &cout );

  tm.toc();
  fmt::print(" {} [ms] done\n",tm.elapsed_ms());
  
  resid = rhs - A*x;

  fmt::print("error    (bicgstab) = {}\n", dist2( x, exact ) );
  fmt::print("residual (bicgstab) = {}\n", normi( resid ) );
#endif

#if 1
  x = 0;
  MA41<double> ma41;
  ma41.load( A );
  fmt::print("\nsolve    (MA41) ... ");
  tm.tic();
  ma41.solve( rhs, x );
  tm.toc();
  fmt::print(" {} [ms] done\n", tm.elapsed_ms());

  resid = rhs - A*x;
  fmt::print("error    (MA41) = {}\n", dist2( x, exact ) );
  fmt::print("residual (MA41) = {}\n", normi( resid ) );
#endif

#if 0
  x = 0;
  MA48<double> ma48;
  ma48.load( A );
  fmt::print("\nsolve    (MA48) ... ");
  tm.tic();
  ma48.solve( rhs, x );
  tm.toc();
  fmt::print(" {} [ms] done\n", tm.elapsed_ms() );

  resid = rhs - A*x;
  fmt::print("error    (MA48) = {}\n", dist2( x, exact ) );
  fmt::print("residual (MA48) = {}\n", normi( resid ) );
#endif

#if 1
  x = 0;
  mkl_PardisoRealU pardiso;
  pardiso.load( A );
  //pardiso.check_matrix();
  pardiso.factorize();
  fmt::print("\nsolve    (pardiso) ... ");
  tm.tic();
  pardiso.solve( rhs, x );
  tm.toc();
  fmt::print(" {} [ms] done\n", tm.elapsed_ms());

  resid = rhs - A*x;
  fmt::print("error    (pardiso) = {}\n", dist2( x, exact ) );
  fmt::print("residual (pardiso) = {}\n", normi( resid ) );
#endif

#if 1
  x = 0;
  SuperLU<double> superlu;
  superlu.load( A );
  fmt::print("\nsolve    (superlu) ... ");
  tm.tic();
  superlu.solve( rhs, x );
  tm.toc();
  fmt::print(" {} [ms] done\n", tm.elapsed_ms());

  resid = rhs - A*x;
  fmt::print("error    (superlu) = {}\n", dist2( x, exact ) );
  fmt::print("residual (superlu) = {}\n", normi( resid ) );
#endif

#if 1
  x = 0;
  UMF<double> umf;
  umf.load( A );
  fmt::print("\nsolve    (UMF) ... ");
  tm.tic();
  umf.solve( rhs, x );
  tm.toc();
  fmt::print(" {} [ms] done\n", tm.elapsed_ms());

  resid = rhs - A*x;
  fmt::print("error    (UMF) = {}\n", dist2( x, exact ) );
  fmt::print("residual (UMF) = {}\n", normi( resid ) );
#endif

  resid = rhs - A*exact;
  fmt::print( "\n\nresidual (exact) = {}\n\n", normi( resid ) );
}

int
main() {
  char const * rMatrix[] = {
    "af23560.mtx.gz", // MA48 fails
    //"memplus.mtx.gz", //
    //"fidap005.mtx.gz",
    //"ASIC_100k.mtx.gz",           // 99340 (ok) 200 iter
    //"ASIC_320ks.mtx.gz",          // 321671 (ok) 200 iter
    //"ASIC_680k.mtx.gz",           // 682862 (ok) 200 iter
    //"af23560.mtx.gz",             // 23560 (no)
    //"audikw_1.mtx.gz",            // 943695 (ok) 11 iter
    //"barrier2-2.mtx.gz",          // 113076 (no)
    //"bmw3_2.mtx.gz",              // 227362 (ok) 8 iter
    //"cage15.mtx.gz",              // In reading Matrix Market File, bad pattern index on line 38025901
    //"CO.mtx.gz",                  // 221119 (ok) 10 iter
    //"dwg961b.mtx.gz",             // 961 (ok) 14
    //"ecology2.mtx.gz",            // 999999 (ok) 53 iter
    //"fidapm05.mtx.gz",            // 42 (no)
    //"memchip.mtx.gz",             // 2707524 (ok) 200 iter 
    //"hor__131.mtx.gz",            // 434 (ok) 200 iter
    //"ldoor.mtx.gz",               // 952203 (ok) 11 iter 
    //"para-9.mtx.gz",              // 155924 (no)
    //"parabolic_fem.mtx.gz",       // 525825 (ok) 54 iter
    //"plat1919.mtx.gz",            // 1919 (ok) 200 iter
    //"s3dkq4m2.mtx.gz",            // 90449 (ok) 10 iter
    nullptr
  };

  for ( char const **p = rMatrix; *p != nullptr; ++p ) {
    string fname = string("mm/")+*p;
    ifstream file( fname.c_str() );
    if ( !file.good() ) {
      cerr << "Cannot open file: " << fname << "\n";
      exit(0);
    }
    igzstream gz(file);
    //ifstream file( fname . c_str() );
    //cout << (file . good()?"OK":"NO") << endl;
    //cout << (file . fail()?"OK":"NO") << endl; 
    testSparseTool( gz );
    file.close();
  }
  fmt::print("\nAll Done Folks\n\n");
  return 0;
}
