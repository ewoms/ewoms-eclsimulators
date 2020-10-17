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
#ifndef EWOMS_TIMESTEPCONTROLINTERFACE_HH
#define EWOMS_TIMESTEPCONTROLINTERFACE_HH

namespace Ewoms
{

    ///////////////////////////////////////////////////////////////////
    ///
    ///  RelativeChangeInterface
    ///
    ///////////////////////////////////////////////////////////////////
    class RelativeChangeInterface
    {
    protected:
        RelativeChangeInterface() {}
    public:
        /// \return || u^n+1 - u^n || / || u^n+1 ||
        virtual double relativeChange() const = 0;

        /// virtual destructor (empty)
        virtual ~RelativeChangeInterface() {}
    };

    ///////////////////////////////////////////////////////////////////
    ///
    ///  TimeStepControlInterface
    ///
    ///////////////////////////////////////////////////////////////////
    class TimeStepControlInterface
    {
    protected:
        TimeStepControlInterface() {}
    public:
        /// compute new time step size suggestions based on the PID controller
        /// \param dt          time step size used in the current step
        /// \param iterations  number of iterations used (linear/nonlinear)
        /// \param timeError   object to compute || u^n+1 - u^n || / || u^n+1 ||
        ///
        /// \return suggested time step size for the next step
        virtual double computeTimeStepSize( const double dt, const int iterations, const RelativeChangeInterface& relativeChange , const double simulationTimeElapsed) const = 0;

        /// virtual destructor (empty)
        virtual ~TimeStepControlInterface () {}
    };

}
#endif
