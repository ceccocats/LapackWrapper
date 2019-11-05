/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  SparseTool   : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH UMFPACK  |
 |                                                                          |
 |  date         : 2008, 14 April                                           |
 |  version      : 1.0                                                      |
 |  file         : test.UMF.cc                                              |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Dipartimento di Ingegneria Industriale                   |
 |                 University of Trento                                     |
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

//#include <sparse_tool/interfaces/MA41.hh>
#include <sparse_tool/interfaces/mkl_pardiso.hh>
#include <sparse_tool/interfaces/Pardiso.hh>

#include <lapack_wrapper/TicToc.hh>

#include <fstream>
#include <iostream>

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshorten-64-to-32"
#endif

#include "zstream/izstream.hh"

#ifdef __clang__
#pragma clang diagnostic pop
#endif

using namespace SparseToolLoad;
using namespace std;
using namespace zstream;

typedef std::complex<double> cplx;

void
testSparseTool( istream & mm_file ) {
  TicToc                   tm;
  MatrixMarket             mm;
  ILDUpreconditioner<cplx> ildu;
  //ILDUKpreconditioner<cplx> ildu;
  CCoorMatrix<cplx>        A;
  Vector<cplx>             x, rhs, exact, resid;

  fmt::print("read matrix...");
  tm.tic();
  mm.read( mm_file, A );
  tm.toc();
  fmt::print(" {} [ms] done\n", tm.elapsed_ms());

  //Spy( mm_file + ".eps" , A, 15.0 );
  
  exact . resize( A.numRows() );
  x     . resize( A.numRows() );
  rhs   . resize( A.numRows() );
  resid . resize( A.numRows() );

  fmt::print("factorize (ildu) ...");
  tm.tic();
  ildu.build(A);
  tm.toc();
  fmt::print(" {} [ms] done\n", tm.elapsed_ms());

  exact = cplx(1);
  rhs   = A * exact;

  fmt::print("solve (ildu) ... ");
  tm.tic();

  double   epsi    = 1e-15;
  unsigned maxIter = 200;
  //unsigned maxSubIter = 50;
  unsigned iter;
  double   res = bicgstab( A, rhs, x, ildu, epsi, maxIter, iter, &cout );
  //double   res = gmres( A, rhs, x, ildu, epsi, maxSubIter, maxIter, iter, &cout );

  tm.toc();
  fmt::print( " {} [s] done\n", tm.elapsed_s());
  
  resid = rhs - A*x;

  fmt::print("\nerror    (ildu) = {}\n", dist2( x, exact ) );
  fmt::print("\nresidual (ildu) = {}\n", normi( resid ) );

#if 1

  tm.tic();
  mkl_PardisoComplexU pardiso;
  pardiso.load( A );
  pardiso.factorize();
  pardiso.solve( rhs, x );

  tm.toc();

  resid = rhs - A*x;

  fmt::print("error    (pardiso) = {}\n", dist2( x, exact ) );
  fmt::print("residual (pardiso) = {}\n", normi( resid ) );
  fmt::print("elapsed time {} [s]\n", tm.elapsed_s() );
#endif

}

int
main() {
  char const * rMatrix[] = {
    //"conf5_4-8x8-05.mtx.gz", // no
    "dielFilterV2clx.mtx.gz", // 607232, ok 52 gmres, 11 bcgstab
    //"dielFilterV3clx.mtx.gz", // 420408, ok 53 gmres, ok 11 bicgstab
    //"fem_filter.mtx.gz", // 74062, fail gmres, fail bcgstab
    //"fem_hifreq_circuit.mtx.gz", // 491100, fail gmres, fail bicgstab
    //"mono_500Hz.mtx.gz", // 169410 fail gmres, fail bicgstab
    //"mplate.mtx.gz", // 5962, 7 bicgstab
    //"qc2534.mtx.gz", // 2534, 35 bicgstab
    //"qc324.mtx.gz", // ok 5 bicgstab
    //"RFdevice.mtx.gz", // 74104 fail gmres, fail bicgstab
    //"ted_AB.mtx.gz", // 10605 fail gmres, fail bicgstab
    //"ted_AB_unscaled.mtx.gz",
    //"vfem.mtx.gz", // 93476 fail gmres, fail bicgstab
    //"windscreen_M.mtx.gz", // 22692 ok 200 gmres, 5 bicgstab
    //"windscreen.mtx.gz", // 22692 ok 4 gmres, 3 bicgstab
    nullptr
  };

  for ( char const **p = rMatrix; *p != nullptr; ++p ) {
    string fname = string("mmc/")+*p;
    ifstream file( fname.c_str() );
    if ( !file.good() ) {
      fmt::print("Cannot open file: {}\n", fname);
      exit(0);
    }
    igzstream gz(file);
    //ifstream file( fname.c_str() );
    //cout << (file.good()?"OK":"NO") << endl;
    //cout << (file.fail()?"OK":"NO") << endl;
    testSparseTool( gz );
    file.close();
  }
  fmt::print("\nAll Done Folks\n\n");
  return 0;
}
