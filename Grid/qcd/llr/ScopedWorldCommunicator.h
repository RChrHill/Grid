/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ScopedWorldCommunicator.h

Copyright (C) 2026

Author: Ryan Hill <Ryan.Hill@ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution
directory
*************************************************************************************/
/*  END LEGAL */

#pragma once

#include <Grid/communicator/Communicator.h>

NAMESPACE_BEGIN(Grid);

// HACK: This is used to temporarily swap over communicators held by each
// replica in order to distribute different RNG results to each replica.
class ScopedWorldCommunicator
{
public:
    ScopedWorldCommunicator(const CartesianCommunicator& scope)
        : saved_(CartesianCommunicator::communicator_world)
    {
        CartesianCommunicator::communicator_world = scope.communicator;
    }

    ~ScopedWorldCommunicator()
    {
        CartesianCommunicator::communicator_world = saved_;
    }

    // Delete copy ops to prevent dtor calls
    ScopedWorldCommunicator(const ScopedWorldCommunicator&) = delete;
    ScopedWorldCommunicator& operator=(const ScopedWorldCommunicator&) = delete;
private:
    Grid_MPI_Comm saved_;
};

NAMESPACE_END(Grid);
