// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
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
#ifndef EWOMS_PRECONDITIONERWITHUPDATE_HH
#define EWOMS_PRECONDITIONERWITHUPDATE_HH

#include <dune/istl/preconditioner.hh>
#include <dune/istl/solvercategory.hh>
#include <memory>
#include <boost/property_tree/ptree.hpp>
namespace Dune
{

/// Interface class adding the update() method to the preconditioner interface.
template <class X, class Y>
class PreconditionerWithUpdate : public Preconditioner<X, Y>
{
public:
    virtual void update() = 0;
};

template <class OriginalPreconditioner>
class DummyUpdatePreconditioner : public PreconditionerWithUpdate<typename OriginalPreconditioner::domain_type,
                                                                  typename OriginalPreconditioner::range_type>
{
public:
    template <class... Args>
    DummyUpdatePreconditioner(Args&&... args)
        : orig_precond_(std::forward<Args>(args)...)
    {
    }

    using X = typename OriginalPreconditioner::domain_type;
    using Y = typename OriginalPreconditioner::range_type;

    virtual void pre(X& x, Y& b) override
    {
        orig_precond_.pre(x, b);
    }

    virtual void apply(X& v, const Y& d) override
    {
        orig_precond_.apply(v, d);
    }

    virtual void post(X& x) override
    {
        orig_precond_.post(x);
    }

    virtual SolverCategory::Category category() const
    {
        return orig_precond_.category();
    }

    // The update() function does nothing for a wrapped preconditioner.
    virtual void update() override
    {
    }

private:
    OriginalPreconditioner orig_precond_;
};

template <class OriginalPreconditioner, class... Args>
std::shared_ptr<DummyUpdatePreconditioner<OriginalPreconditioner>>
wrapPreconditioner(Args&&... args)
{
    return std::make_shared<DummyUpdatePreconditioner<OriginalPreconditioner>>(std::forward<Args>(args)...);
}

} // namespace Dune

#endif // EWOMS_PRECONDITIONERWITHUPDATE_HH
