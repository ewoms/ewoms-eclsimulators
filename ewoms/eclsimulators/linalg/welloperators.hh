/*

  This file is part of the eWoms project.

  eWoms is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  eWoms is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with eWoms.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef EWOMS_WELLOPERATORS_HH
#define EWOMS_WELLOPERATORS_HH

#include <dune/istl/operators.hh>

namespace Ewoms
{

//=====================================================================
// Implementation for ISTL-matrix based operators
// Note: the classes WellModelMatrixAdapter and
// WellModelGhostLastMatrixAdapter were moved from ISTLSolver.hpp
// and subsequently modified.
//=====================================================================

/// Linear operator wrapper for well model.
///
/// This class is intended to hide the actual type of the well model
/// (which depends on a TypeTag) by making a simple linear operator
/// wrapper. That way the WellModelMatrixAdapter template does not need
/// the concrete WellModel type, and we can avoid instantiating
/// WellModelMatrixAdapter separately for each TypeTag, as it will only
/// depend on the matrix and vector types involved, which typically are
/// just one for each block size with block sizes 1-4.
template <class WellModel, class X, class Y>
class WellModelAsLinearOperator : public Dune::LinearOperator<X, Y>
{
public:
    using Base = Dune::LinearOperator<X, Y>;
    using field_type = typename Base::field_type;

    explicit WellModelAsLinearOperator(const WellModel& wm)
        : wellMod_(wm)
    {
    }
    /*! \brief apply operator to x:  \f$ y = A(x) \f$
          The input vector is consistent and the output must also be
       consistent on the interior+border partition.
     */
    void apply (const X& x, Y& y) const override
    {
        wellMod_.apply(x, y );
    }

    //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
    virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const override
    {
        wellMod_.applyScaleAdd( alpha, x, y );
    }

    /// Category for operator.
    /// This is somewhat tricky, I consider this operator sequential
    /// since (unlike WellModelMatrixAdapter) it does not do any
    /// projections etc. but only forwards the calls to the sequential
    /// well model.
    Dune::SolverCategory::Category category() const override
    {
        return Dune::SolverCategory::sequential;
    }
private:
    const WellModel& wellMod_;
};

/*!
   \brief Adapter to combine a matrix and another linear operator into
   a combined linear operator.

   Adapts a matrix A plus another linear operator W (typically from
   wells) to the assembled linear operator interface by returning S
   from getmat() and making apply() and applyscaleadd() apply both A
   and W to the input vector. In addition this is a parallel-aware
   adapter, that does not require the W operator to be parallel, but
   makes it into one by making the proper projections.
 */
template<class M, class X, class Y, bool overlapping >
class WellModelMatrixAdapter : public Dune::AssembledLinearOperator<M,X,Y>
{
public:
  typedef M matrix_type;
  typedef X domain_type;
  typedef Y range_type;
  typedef typename X::field_type field_type;

#if HAVE_MPI
  typedef Dune::OwnerOverlapCopyCommunication<int,int> communication_type;
#else
  typedef Dune::CollectiveCommunication< int > communication_type;
#endif

  Dune::SolverCategory::Category category() const override
  {
    return overlapping ?
           Dune::SolverCategory::overlapping : Dune::SolverCategory::sequential;
  }

  //! constructor: just store a reference to a matrix
  WellModelMatrixAdapter (const M& A,
                          const Dune::LinearOperator<X, Y>& wellOper,
                          const std::shared_ptr< communication_type >& comm = std::shared_ptr< communication_type >())
      : A_( A ), wellOper_( wellOper ), comm_(comm)
  {}

  virtual void apply( const X& x, Y& y ) const override
  {
    A_.mv( x, y );

    // add well model modification to y
    wellOper_.apply(x, y );

#if HAVE_MPI
    if( comm_ )
      comm_->project( y );
#endif
  }

  // y += \alpha * A * x
  virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const override
  {
    A_.usmv(alpha,x,y);

    // add scaled well model modification to y
    wellOper_.applyscaleadd( alpha, x, y );

#if HAVE_MPI
    if( comm_ )
      comm_->project( y );
#endif
  }

  virtual const matrix_type& getmat() const override { return A_; }

protected:
  const matrix_type& A_ ;
  const Dune::LinearOperator<X, Y>& wellOper_;
  std::shared_ptr< communication_type > comm_;
};

/*!
   \brief Adapter to combine a matrix and another linear operator into
   a combined linear operator.

   This is similar to WellModelMatrixAdapter, with the difference that
   here we assume a parallel ordering of rows, where ghost rows are
   located after interior rows.
 */
template<class M, class X, class Y, bool overlapping >
class WellModelGhostLastMatrixAdapter : public Dune::AssembledLinearOperator<M,X,Y>
{
public:
    typedef M matrix_type;
    typedef X domain_type;
    typedef Y range_type;
    typedef typename X::field_type field_type;

#if HAVE_MPI
    typedef Dune::OwnerOverlapCopyCommunication<int,int> communication_type;
#else
    typedef Dune::CollectiveCommunication< int > communication_type;
#endif

    Dune::SolverCategory::Category category() const override
    {
        return overlapping ?
            Dune::SolverCategory::overlapping : Dune::SolverCategory::sequential;
    }

    //! constructor: just store a reference to a matrix
    WellModelGhostLastMatrixAdapter (const M& A,
                                     const Dune::LinearOperator<X, Y>& wellOper,
                                     const size_t interiorSize )
        : A_( A ), wellOper_( wellOper ), interiorSize_(interiorSize)
    {}

    virtual void apply( const X& x, Y& y ) const override
    {
        for (auto row = A_.begin(); row.index() < interiorSize_; ++row)
        {
            y[row.index()]=0;
            auto endc = (*row).end();
            for (auto col = (*row).begin(); col != endc; ++col)
                (*col).umv(x[col.index()], y[row.index()]);
        }

        // add well model modification to y
        wellOper_.apply(x, y );

        ghostLastProject( y );
    }

    // y += \alpha * A * x
    virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const override
    {
        for (auto row = A_.begin(); row.index() < interiorSize_; ++row)
        {
            auto endc = (*row).end();
            for (auto col = (*row).begin(); col != endc; ++col)
                (*col).usmv(alpha, x[col.index()], y[row.index()]);
        }
        // add scaled well model modification to y
        wellOper_.applyscaleadd( alpha, x, y );

        ghostLastProject( y );
    }

    virtual const matrix_type& getmat() const override { return A_; }

protected:
    void ghostLastProject(Y& y) const
    {
        size_t end = y.size();
        for (size_t i = interiorSize_; i < end; ++i)
            y[i] = 0;
    }

    const matrix_type& A_ ;
    const Dune::LinearOperator<X, Y>& wellOper_;
    size_t interiorSize_;
};

} // namespace Ewoms

#endif // EWOMS_WELLOPERATORS_HH

