/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  SparseTool   : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH SuperLU  |
 |                                                                          |
 |  date         : 2008, 14 April                                           |
 |  version      : 1.0                                                      |
 |  file         : SparseTool_SuperLU.hh                                    |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Dipartimento di Ingegneria Industriale                   |
 |                 Universita` degli Studi di Trento                        |
 |                 email: enrico.bertolazzi@unitn.it                        |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifndef SPARSETOOL_SUPERLU_HH
#define SPARSETOOL_SUPERLU_HH

#include "SparseTool.hh"

#include <complex>

/* 
 * SUPERLU interface http://crd.lbl.gov/~xiaoye/SuperLU/
 */

namespace SparseTool {

  template <typename T> struct RealType { typedef T Type; };

  template <> struct RealType<std::complex<float> >  { typedef float Type; };
  template <> struct RealType<std::complex<double> > { typedef double Type; };

  template <typename T>
  class superLU {
  public:
    typedef T valueType;
    typedef typename RealType<T>::Type realType;

  private:
    void * L;
    void * U;
    int  * perm_c; // row permutations from partial pivoting
    int  * perm_r; // column permutation vector

    bool     allocated;
    unsigned Lnnz;
    unsigned Unnz;
    unsigned mem_usage_for_lu;
    unsigned mem_usage_total;

    void freeMemory();

  public:

    superLU() : allocated(false) {};
    ~superLU() { freeMemory(); }

    int
    load( CColMatrix<T> const & A, realType const dropTolerance = 0 );

    template <typename MT>
    int
    load( Sparse<T,MT> const & A, realType const dropTolerance = 0 ) {
      CColMatrix<T> A1 = A;
      return load(A1,dropTolerance);
    }

    int
    solve( Vector<T> const & b,
           Vector<T>       & x,
           bool transpose = false );

    void
    statistic( ostream & s ) {
      s <<   "No of nonzeros in factor L = " << Lnnz
        << "\nNo of nonzeros in factor U = " << Unnz
        << "\nNo of nonzeros in L+U      = " << Lnnz + Unnz
        << "\nMemory usage for LU (MB)   = " << mem_usage_for_lu/1E6  
        << "\nTotal memory usage  (MB)   = " << mem_usage_total/1E6
        << "\n";
    }

  };

  /*
  //  ###  #       #     #
  //   #   #       #     #
  //   #   #       #     #
  //   #   #       #     #
  //   #   #       #     #
  //   #   #       #     #
  //  ###  #######  #####
  */
  /*! \class SLUPreco
      \brief Incomplete \c LU preconditioner
   */
  template <typename T>
  class superLUpreconditioner : public Preco<superLUpreconditioner<T> > {
  public:

    // \cond NODOC
    typedef superLUpreconditioner<T> SLUPRECO;
    typedef Preco<SLUPRECO>          PRECO;
    // \endcond
    typedef T valueType; //!< type of the element of the preconditioner
    typedef typename RealType<T>::Type realType;

  private:

    mutable superLU<T> ILU;

  public:

    superLUpreconditioner(void) : Preco<SLUPRECO>() {}
    
    template <typename MAT>
    superLUpreconditioner( MAT const & M, realType const dropTolerance ) : Preco<SLUPRECO>() 
    { ILU.load( M, dropTolerance ); }

    //! build the preconditioner from matrix \c M
    template <typename MAT>
    void
    build( MAT const & M, realType const dropTolerance )
    { ILU.load( M, dropTolerance ); }

    //! apply preconditioner to vector \c v and store result to vector \c res
    template <typename VECTOR>
    void
    assPreco( VECTOR & res, VECTOR const & v ) const
    { ILU.solve( v, res, false ); }

  };

  //! \cond NODOC
  template <typename T, typename TP> inline
  Vector_V_div_P<Vector<T>,superLUpreconditioner<TP> >
  operator / (Vector<T> const & v, superLUpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,superLUpreconditioner<TP> >(v,P); }
  //! \endcond

}

namespace SparseToolLoad {
  using ::SparseTool::superLU;
  using ::SparseTool::superLUpreconditioner;
}

#endif
