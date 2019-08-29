/*--------------------------------------------------------------------------*\
 |                                                                          |
 |  SparseTool   : DRIVER FOR TESTING THE TOOLKIT INTERFACING WITH UMFPACK  |
 |                                                                          |
 |  date         : 2008, 6 May                                              |
 |  version      : 1.0.1                                                    |
 |  file         : SparseTool_Umf.hh                                        |
 |  authors      : Enrico Bertolazzi                                        |
 |  affiliations : Dipartimento di Ingegneria Industriale                   |
 |                 Universita` degli Studi di Trento                        |
 |                 email: enrico.bertolazzi@unitn.it                        |
 |                                                                          |
\*--------------------------------------------------------------------------*/

#ifndef SPARSETOOL_UMF_HH
#define SPARSETOOL_UMF_HH

#include "SparseTool.hh"
#include <umfpack.h>

#include <complex>

namespace SparseTool {

  template <typename T>
  class UMF {
  public:  
    typedef T valueType;
    typedef typename return_trait<T>::valueType realType;

  private:

    realType Info    [UMFPACK_INFO];
    realType Control [UMFPACK_CONTROL];

    CColMatrix<T>    A_buffer;

    indexType                 numRow;
    Vector<valueType> const * Ax;
    Vector<indexType> const * Ai;
    Vector<indexType> const * Ap;
    
    Vector<realType> Areal, Aimag;

    void *Symbolic, *Numeric;
    bool allocated;

    int  load( realType const dropTolance );
    void free();

  public:
  
    UMF() : allocated(false) {};
    ~UMF() { free(); }

    int load( CCoorMatrix<T> const & Mat, realType const dropTolerance = 0 );
    int load( CColMatrix<T>  const & Mat, realType const dropTolerance = 0 );

    template <typename MT>
    int
    load( Sparse<T,MT> const & A, realType const dropTolerance = 0 ) {
      A_buffer = A;
      return load(A_buffer,dropTolerance);
    }

    int
    solve( Vector<T> const & b,
           Vector<T>       & x,
           bool transpose = false );

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
  /*! \class UMFpreconditioner
      \brief Incomplete \c LU preconditioner
   */
  template <typename T>
  class UMFpreconditioner : public Preco<UMFpreconditioner<T> > {
  public:

    // \cond NODOC
    typedef UMFpreconditioner<T> UMFPRECO;
    typedef Preco<UMFPRECO>      PRECO;
    // \endcond
    typedef T valueType; //!< type of the element of the preconditioner
    typedef typename return_trait<T>::valueType realType;

  private:

    mutable UMF<T> ILU;

  public:

    UMFpreconditioner(void) : Preco<UMFPRECO>() {}
    
    template <typename MAT>
    UMFpreconditioner( MAT const & M, realType const dropTolerance ) : Preco<UMFPRECO>() 
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
  Vector_V_div_P<Vector<T>,UMFpreconditioner<TP> >
  operator / (Vector<T> const & v, UMFpreconditioner<TP> const & P)
  { return Vector_V_div_P<Vector<T>,UMFpreconditioner<TP> >(v,P); }
  //! \endcond

}

namespace SparseToolLoad {
  using ::SparseTool::UMF;
  using ::SparseTool::UMFpreconditioner;
}

#endif
