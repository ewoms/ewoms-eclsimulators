/*
  The content of this file is based on the file dune/istl/matrixmarket.hh in
  the Dune module dune-istl.

  The license of this file is therefore the same as that of Dune, see
  https://www.dune-project.org/about/license/
*/

#ifndef EWOMS_MATRIXMARKETSPECIALIZATIONS_HH
#define EWOMS_MATRIXMARKETSPECIALIZATIONS_HH

#include <dune/istl/matrixmarket.hh>

namespace Ewoms
{
template<typename T, int i, int j>
class MatrixBlock;
}

namespace Dune
{

namespace MatrixMarketImpl
{

    template <typename T, int i, int j, typename A>
    struct mm_header_printer<BCRSMatrix<Ewoms::MatrixBlock<T,i,j>, A>>
    {
        static void print(std::ostream& os)
        {
            os << "%%MatrixMarket matrix coordinate ";
            os << mm_numeric_type<T>::str() << " general" << std::endl;
        }
    };

    template <typename T, int i, int j, typename A>
    struct mm_block_structure_header<BCRSMatrix<Ewoms::MatrixBlock<T,i,j>, A>>
    {
        using M = BCRSMatrix<Ewoms::MatrixBlock<T,i,j>, A>;
        static void print(std::ostream& os, const M&)
        {
            os << "% ISTL_STRUCT blocked ";
            os << i << " " << j << std::endl;
        }
    };

} // namespace MatrixMarketImpl

template<class M>
struct mm_multipliers;

template <typename T, int i, int j, typename A>
struct mm_multipliers<BCRSMatrix<Ewoms::MatrixBlock<T,i,j>, A>>
{
    enum {
        rows = i,
        cols = j
    };
};

} // namespace Dune

#endif
