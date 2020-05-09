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
#ifndef EWOMS_PARALLELRESTRICTEDADDITIVESCHWARZ_HH
#define EWOMS_PARALLELRESTRICTEDADDITIVESCHWARZ_HH

#include <dune/istl/preconditioner.hh>
#include <dune/istl/paamg/smoother.hh>

namespace Ewoms
{

template<class X, class Y, class C, class T>
class ParallelRestrictedOverlappingSchwarz;

} // end namespace Ewoms

namespace Dune
{

namespace Amg
{

/// \brief Tells AMG how to construct the Ewoms::ParallelOverlappingILU0 smoother
/// \tparam Domain The type of the Vector representing the domain.
/// \tparam Range The type of the Vector representing the range.
/// \tparam ParallelInfo The type of the parallel information object
///         used, e.g. Dune::OwnerOverlapCommunication
/// \tparam SeqPreconditioner The underlying sequential preconditioner to use.
template<class Range, class Domain, class ParallelInfo, class SeqPreconditioner>
class ConstructionTraits<Ewoms::ParallelRestrictedOverlappingSchwarz<Range,
                                                                    Domain,
                                                                    ParallelInfo,
                                                                    SeqPreconditioner> >
{
public:
    typedef DefaultParallelConstructionArgs<SeqPreconditioner,ParallelInfo> Arguments;
    typedef ConstructionTraits<SeqPreconditioner> SeqConstructionTraits;

    /// \brief Construct a parallel restricted overlapping schwarz preconditioner.
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
    typedef std::shared_ptr< Ewoms::ParallelRestrictedOverlappingSchwarz<Range,
                                                                       Domain,
                                                                       ParallelInfo,
                                                                       SeqPreconditioner> > ParallelRestrictedOverlappingSchwarzPointer;
#else
    typedef Ewoms::ParallelRestrictedOverlappingSchwarz<Range,
                                                      Domain,
                                                      ParallelInfo,
                                                      SeqPreconditioner>*  ParallelRestrictedOverlappingSchwarzPointer;
#endif

    static inline ParallelRestrictedOverlappingSchwarzPointer
    construct(Arguments& args)
    {
        return
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
             std::make_shared(
#endif
                new Ewoms::ParallelRestrictedOverlappingSchwarz
                              <Range,Domain,ParallelInfo,SeqPreconditioner>(*SeqConstructionTraits ::construct(args),
                                                         args.getComm())
#if DUNE_VERSION_NEWER(DUNE_ISTL, 2, 7)
        );
#else
        ;
#endif
    }

    /// \brief Deconstruct and free a parallel restricted overlapping schwarz preconditioner.
    static inline void deconstruct(Ewoms::ParallelRestrictedOverlappingSchwarz
                                   <Range,Domain,ParallelInfo,SeqPreconditioner>* bp)
    {
        SeqConstructionTraits
            ::deconstruct(static_cast<SeqPreconditioner*>(&bp->preconditioner));
        delete bp;
    }

};

/// \brief Tells AMG how to use Ewoms::ParallelOverlappingILU0 smoother
/// \tparam Domain The type of the Vector representing the domain.
/// \tparam Range The type of the Vector representing the range.
/// \tparam ParallelInfo The type of the parallel information object
///         used, e.g. Dune::OwnerOverlapCommunication
/// \tparam SeqPreconditioner The underlying sequential preconditioner to use.
template<class Range, class Domain, class ParallelInfo, class SeqPreconditioner>
struct SmootherTraits<Ewoms::ParallelRestrictedOverlappingSchwarz<Range,
                                                                Domain,
                                                                ParallelInfo,
                                                                SeqPreconditioner> >
{
    typedef DefaultSmootherArgs<typename SeqPreconditioner::matrix_type::field_type> Arguments;

};

} // end namespace Amg

} // end namespace Dune

namespace Ewoms{

/// \brief Block parallel preconditioner.
///
/// This is essentially a wrapper that takes a sequential
/// preconditioner. In each step the sequential preconditioner
/// is applied to the whole subdomain and then all owner data
/// points are updated on all other processes from the processor
/// that knows the complete matrix row for this data point (in dune-istl
/// speak that is the one that owns the data).
///
/// Note that this is different from the usual approach in dune-istl where
/// the application of the sequential preconditioner only takes place on
/// the (owner) partition of the process disregarding any overlap/ghost region.
///
/// For more information see https://www.cs.colorado.edu/~cai/papers/rash.pdf
///
/// \tparam Domain The type of the Vector representing the domain.
/// \tparam Range The type of the Vector representing the range.
/// \tparam ParallelInfo The type of the parallel information object
///         used, e.g. Dune::OwnerOverlapCommunication
/// \tparam SeqPreconditioner The underlying sequential preconditioner to use.
template<class Range, class Domain, class ParallelInfo, class SeqPreconditioner=Dune::Preconditioner<Range,Domain> >
class ParallelRestrictedOverlappingSchwarz
    : public Dune::Preconditioner<Range,Domain> {
    friend class Dune::Amg
    ::ConstructionTraits<ParallelRestrictedOverlappingSchwarz<Range,
                                                              Domain,
                                                              ParallelInfo,
                                                              SeqPreconditioner> >;
public:
    //! \brief The domain type of the preconditioner.
    typedef Domain domain_type;
    //! \brief The range type of the preconditioner.
    typedef Range range_type;
    //! \brief The field type of the preconditioner.
    typedef typename Domain::field_type field_type;
    //! \brief The type of the communication object.
    typedef ParallelInfo communication_type;

    // define the category
    enum {
        //! \brief The category the precondtioner is part of.
        category=Dune::SolverCategory::overlapping
    };

    /*! \brief Constructor.

      constructor gets all parameters to operate the prec.
      \param p The sequential preconditioner.
      \param c The communication object for syncing overlap and copy
      data points. (E.~g. OwnerOverlapCommunication )
    */
    ParallelRestrictedOverlappingSchwarz (SeqPreconditioner& p, const communication_type& c)
        : preconditioner_(p), communication_(c)
    {   }

    /*!
      \brief Prepare the preconditioner.

      \copydoc Preconditioner::pre(X&,Y&)
    */
    virtual void pre (Domain& x, Range& b)
    {
        communication_.copyOwnerToAll(x,x);     // make dirichlet values consistent
        preconditioner_.pre(x,b);
    }

    /*!
      \brief Apply the preconditioner

      \copydoc Preconditioner::apply(X&,const Y&)
    */
    virtual void apply (Domain& v, const Range& d)
    {
        apply<true>(v, d);
    }

    template<bool forward>
    void apply (Domain& v, const Range& d)
    {
        // hack us a mutable d to prevent copying.
        Range& md = const_cast<Range&>(d);
        communication_.copyOwnerToAll(md,md);
        preconditioner_.template apply<forward>(v,d);
        communication_.copyOwnerToAll(v,v);
        // Make sure that d is the same as at the beginning of apply.
        communication_.project(md);
    }

    /*!
      \brief Clean up.

      \copydoc Preconditioner::post(X&)
    */
    virtual void post (Range& x)
    {
        preconditioner_.post(x);
    }

private:
    //! \brief a sequential preconditioner
    SeqPreconditioner& preconditioner_;

    //! \brief the communication object
    const communication_type& communication_;
};

} // end namespace Ewoms
#endif
