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
#ifndef EWOMS_OWNINGBLOCKPRECONDITIONER_HH
#define EWOMS_OWNINGBLOCKPRECONDITIONER_HH

#include <ewoms/eclsimulators/linalg/preconditionerwithupdate.hh>

#include <dune/istl/schwarz.hh>

namespace Dune
{

template <class OriginalPreconditioner, class Comm>
class OwningBlockPreconditioner : public PreconditionerWithUpdate<typename OriginalPreconditioner::domain_type,
                                                                  typename OriginalPreconditioner::range_type>
{
public:
    template <class... Args>
    OwningBlockPreconditioner(const Comm& comm, Args&&... args)
        : orig_precond_(std::forward<Args>(args)...)
        , block_precond_(orig_precond_, comm)
    {
    }

    using X = typename OriginalPreconditioner::domain_type;
    using Y = typename OriginalPreconditioner::range_type;

    virtual void pre(X& x, Y& b) override
    {
        block_precond_.pre(x, b);
    }

    virtual void apply(X& v, const Y& d) override
    {
        block_precond_.apply(v, d);
    }

    virtual void post(X& x) override
    {
        block_precond_.post(x);
    }

    virtual SolverCategory::Category category() const
    {
        return block_precond_.category();
    }

    // The update() function does nothing for a wrapped preconditioner.
    virtual void update() override
    {
        orig_precond_.update();
    }

private:
    OriginalPreconditioner orig_precond_;
    BlockPreconditioner<X, Y, Comm, OriginalPreconditioner> block_precond_;
};

template <class OriginalPreconditioner, class Comm, class... Args>
std::shared_ptr<OwningBlockPreconditioner<OriginalPreconditioner, Comm>>
wrapBlockPreconditioner(const Comm& comm, Args&&... args)
{
    return std::make_shared<OwningBlockPreconditioner<OriginalPreconditioner, Comm>>(comm, std::forward<Args>(args)...);
}

} // namespace Dune

#endif // EWOMS_OWNINGBLOCKPRECONDITIONER_HH
