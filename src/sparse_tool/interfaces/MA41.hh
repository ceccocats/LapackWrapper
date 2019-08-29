/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  SparseTool   : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH MA41     |
 |                                                                          |
 |  date         : 2011, 18 July                                            |
 |  version      : 1.0.                                                     |
 |  file         : MA41.hxx                                                 |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Dipartimento di Ingegneria Industriale                   |
 |                 Universita` degli Studi di Trento                        |
 |                 email: enrico.bertolazzi@unitn.it                        |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifndef SPARSETOOL_MA41_HH
#define SPARSETOOL_MA41_HH

#include "sparse_tool.hh"
#include "../../HSL/hsl.h"

namespace SparseTool {

  template <typename T>
  class MA41 {
  public:
    typedef T            real_type;
    typedef HSL::integer integer;

  private:
    enum {
      ANALISYS                             = 1,
      FACTORIZATION                        = 2,
      SOLVE                                = 3,
      ANALISYS_and_FACTORIZATION           = 4,
      FACTORIZATION_and_SOLVE              = 5,
      ANALISYS_and_FACTORIZATION_and_SOLVE = 6
    };

    bool      verbose;
    integer   N, NE;
    integer   KEEP[50];
    integer   ICNTL[20];
    integer   INFO[20];
    real_type CNTL[10];
    real_type RINFO[20];

    // work arrays
    vector<integer>   IS;
    vector<real_type> COLSCA, ROWSCA, S;

    // the matrix
    vector<real_type> A;
    vector<integer>   IRN, JCN;

    void
    msg_info() const {
      std::cout
        << "\nReal    space needs for factors    = " << INFO[2]
        << "\nInteger space needs for factors    = " << INFO[3]
        << "\nMaximum frontal size               = " << INFO[4]
        << "\nNodes in the tree                  = " << INFO[5]
        << "\nMinimum MAXIS                      = " << INFO[6]
        << "\nMinimum MAXS                       = " << INFO[7]
        << "\nReal    Space for LU factorizarion = " << INFO[8]
        << "\nInteger Space for LU factorizarion = " << INFO[9]
        << "\nInteger Space for LU factorizarion = " << INFO[10]
        << "\nLargest frontal matrix             = " << INFO[11]
        << "\n";
    }

    void
    msg_infor() const {
      std::cout
        << "\nNFLOPs for LU factorization    = " << RINFO[0]
        << "\nNFLOPs for assembly process    = " << RINFO[1]
        << "\nNFLOPs for elimination process = " << RINFO[2];
      if ( ICNTL[10] > 0 )
        std::cout
          << "\nInfinity norm of the matrix    = " << RINFO[3]
          << "\nInfinity norm of the solution  = " << RINFO[4]
          << "\nNorm of scaled residual        = " << RINFO[5];
      std::cout << "\n";
    }

    void
    msg_error() const {
      if ( INFO[0] >= 0 ) return;
      std::cerr
        << "\n***** ERROR *****"
        << "\ndim = " << N
        << "\nnnz = " << NE
        << "\n";
      switch ( INFO[0] ) {
      case -1:
        std::cerr << "Value of N is out of range N = " << INFO[1] << "\n";
        break;
      case -2:
        std::cerr << "Value of NE is out of range NE = " << INFO[1] << "\n";
        break;
      case -3:
        std::cerr << "JOB has wrong value or analisys was not performed prior to factorization. JOB = " << INFO[1] << "\n";
        break;
      case -4:
        std::cerr << "Error in permutation error\n";
        break;
      case -5:
        std::cerr << "Not enought space to preprocess the input matrix\n"
                  << "MAXS should be increased to at least " << INFO[1] << "\n";
        break;
      case -6:
        std::cerr << "The matrix is structurally singular\n"
                  << "Estimated rank = " << INFO[1] << "\n";
        break;
      case -7:
        std::cerr << "Error from analysis\n"
                  << "MAXIS should be increased to at least " << INFO[1] << "\n";
        break;
      case -8:
        std::cerr << "Error from numerical factorization\n"
                  << "MAXIS should be increased to at least " << INFO[1] << "\n";
        break;
      case -9:
        std::cerr << "Error from numerical factorization\n"
                  << "MAXS should be increased to at least " << INFO[1] << "\n";
        break;
      case -10:
        std::cerr << "The matrix is numerically singular\n"
                  << "Estimated rank = " << INFO[1] << "\n";
        break;
      case -11:
        std::cerr << "Error from the solution phase\n"
                  << "MAXS should be increased to at least " << INFO[1] << "\n";
        break;
      case -12:
        std::cerr << "Not enought space to postprocess the solution\n"
                  << "MAXS should be increased to at least " << INFO[1] << "\n";
        break;
      }
      std::cerr << "***** ERROR *****\n";
    }

    // private functionalities
    void
    load_matrix( integer const nr ) {
      if ( verbose )
        std::cout << "ma41::load_matrix(...)\n";
      N  = nr;
      NE = integer(A.size());

      COLSCA . resize(N);
      ROWSCA . resize(N);
      IS     . resize(2 * NE + 11 * N + 1);
      S      . resize(2 * NE + 11 * N + 1);

      fill(COLSCA . begin(), COLSCA . end(), 0 );
      fill(ROWSCA . begin(), ROWSCA . end(), 0 );
      fill(IS     . begin(), IS     . end(), 0 );
      fill(S      . begin(), S      . end(), 0 );

      // set pars
      HSL::ma41i<real_type>(CNTL, ICNTL, KEEP);

      if ( verbose )
        std::cout << "done\n";
    }

    void
    symbfac() {
      if ( verbose )
        std::cout << "ma41::symbfac() perform analisys\n";

      // analysis phase
      vector<integer> tmpIS(2 * NE + 12 * N + 1);
      integer JOB = ANALISYS;
      HSL::ma41a<T>(
        JOB, N, NE, &IRN.front(), &JCN.front(), &A.front(),
        nullptr, &COLSCA.front(), &ROWSCA.front(), KEEP,
        &tmpIS.front(), integer(tmpIS.size()),
        &S.front(), integer(S.size()),
        CNTL, ICNTL, INFO, RINFO
      );

      if ( verbose ) msg_error();

      unsigned mem = 10*INFO[6]*sizeof(integer) + 10*INFO[7]*sizeof(T);
      if ( mem > 4*1E9 ) {
        std::cerr << "out of memory mem = " << mem << '\n';
        exit(1);
      }

      IS . resize(10*INFO[6]);
      S  . resize(10*INFO[7]);

      copy( tmpIS.begin(), tmpIS.end(), IS.begin() );

      if ( verbose )
        std::cout << "ma41_wrapper::symbfac() do factorization\n";

      // factorization phase
      // -------------------
      JOB = FACTORIZATION;
      HSL::ma41a<real_type>(
        JOB, N, NE, &IRN.front(), &JCN.front(), &A.front(),
        nullptr, &COLSCA.front(), &ROWSCA.front(), KEEP,
        &IS.front(), integer(IS.size()), &S.front(), integer(S.size()),
        CNTL, ICNTL, INFO, RINFO
      );

      if ( verbose ) {
        msg_info();
        msg_error();
        std::cout << "done\n";
      }

    }

  public:
  
    MA41() {}
    ~MA41() {}

    void
    init(void) {
      IRN . clear();
      JCN . clear();
      A   . clear();
    }

    void
    insert(
      integer   i,
      integer   j,
      real_type a
    ) {
      SPARSETOOL_ASSERT(
        i >= 0 && i < this->N && j >= 0 && j < this->N,
        "MA41::insert( " << i << ", " << j << ", a ) out of matrix"
      )
      IRN . push_back(i+1);
      JCN . push_back(j+1);
      A   . push_back(a);
    }

    void
    setup( integer nr, bool v = false ) {
      verbose = v;
      load_matrix(nr);
      symbfac();
    }

    void
    solve( real_type RHS[] ) {
      integer JOB = SOLVE;
      HSL::ma41a<real_type>(
        JOB, N, NE, &IRN.front(), &JCN.front(), &A.front(),
        RHS, &COLSCA.front(), &ROWSCA.front(), KEEP,
        &IS.front(), integer(IS.size()), &S.front(), integer(S.size()),
        CNTL, ICNTL, INFO, RINFO
      );

      if ( verbose ) {
        msg_infor();
        msg_error();
      }
    }

    template <typename MT>
    int
    load( Sparse<T,MT> const & Mat ) {
      this->init();
      IRN . reserve( Mat.nnz() );
      JCN . reserve( Mat.nnz() );
      A   . reserve( Mat.nnz() );
      this->N = Mat.numRows();
      for ( Mat.Begin(); Mat.End(); Mat.Next() )
        this->insert( Mat.row(), Mat.column(), Mat.value() );
      this->setup( Mat.numRows() );
      return 0;
    }

    int
    solve(
      Vector<T> const & b,
      Vector<T>       & x,
      bool transpose = false
    ) {
      x = b;
      this->solve( &x.front() );
      return 0;
    }

  };
}

namespace SparseToolLoad {
  using ::SparseTool::MA41;
}

#endif
