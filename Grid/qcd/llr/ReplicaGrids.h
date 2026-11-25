/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ReplicaGrids.h

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

#include <Grid/qcd/hmc/HMCModules.h>
#include <Grid/qcd/utils/SpaceTimeGrid.h>

NAMESPACE_BEGIN(Grid);

// Split a Grid into `replica_count` replicas (subgrids).
inline Coordinate createReplicaExchangeGrid(const Coordinate &processors,
                                            int replica_count)
{
    assert(replica_count > 0);
    Coordinate replica_processors = processors;
    int remaining = replica_count;
    for (int d = processors.size() - 1; d >= 0 && remaining > 1; --d)
    {
        int factor = std::gcd(replica_processors[d], remaining);
        replica_processors[d] /= factor;
        remaining /= factor;
    }
    assert(remaining == 1 && "replica_count must divide the processor grid");
    return replica_processors;
}

// This module is used to register the split grid with the HMC
// to be used by the replica exchange module.
inline GridModule ReplicaGridModule(GridCartesian &full_grid,
                                    int replica_count, int &replica_index)
{
  GridModule module;
  module.set_full(new GridCartesian(
      full_grid.FullDimensions(),
      GridDefaultSimd(full_grid.Nd(), vComplex::Nsimd()),
      createReplicaExchangeGrid(full_grid.ProcessorGrid(), replica_count),
      full_grid, replica_index));
  module.set_rb(SpaceTimeGrid::makeFourDimRedBlackGrid(module.get_full()));
  return module;
}

NAMESPACE_END(Grid);
