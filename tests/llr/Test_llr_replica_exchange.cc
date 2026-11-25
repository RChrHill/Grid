/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Test_llr_replica_exchange.cc

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


// MPI test of the replica exchange.
// Run with e.g. mpirun -np 4 ./Test_llr_replica_exchange --grid 4.4.4.4 --mpi 1.1.2.2

#include <Grid/Grid.h>
#include <Grid/qcd/llr/ReplicaExchangeModule.h>
#include <Grid/qcd/llr/ReplicaGrids.h>

using namespace Grid;

int main(int argc, char **argv)
{
    Grid_init(&argc, &argv);

    Coordinate lattice = GridDefaultLatt();
    Coordinate simd = GridDefaultSimd(lattice.size(), vComplex::Nsimd());
    GridCartesian full_grid(lattice, simd, GridDefaultMpi());

    const int replica_count = 2;
    if (full_grid.ProcessorCount() < 2 * replica_count || full_grid.ProcessorCount() % replica_count != 0)
    {
        std::cout << GridLogMessage << "Skipping: needs an even number of ranks, at least 4" << std::endl;
        Grid_finalize();
        return 0;
    }
    int replica_index = -1;
    GridCartesian split_grid(lattice, simd, ReplicaProcessorGrid(full_grid.ProcessorGrid(), replica_count), full_grid, replica_index);

    GridSerialRNG exchange_rng;
    exchange_rng.SeedFixedIntegers({91, 92, 93, 94});

    {
        ScopedWorldCommunicator replica_scope(split_grid);

        GridParallelRNG parallel_rng(&split_grid);
        parallel_rng.SeedFixedIntegers({11, 22, 33, 44, replica_index});
        GridSerialRNG serial_rng;
        serial_rng.SeedFixedIntegers({55, 66, 77, 88, replica_index});

        // Check 1: replicas seeded differently draw different numbers. This checks the
        // ScopedWorldCommunicator is correctly blocking off the different replicas from
        // exchanging RNG values.
        RealD draw = 0.0;
        random(serial_rng, draw);
        std::vector<RealD> draws(replica_count, 0.0); // Same pattern as in ReplicaExchangeModule.h; c.f.
        if (split_grid.IsBoss())
        {
            draws[replica_index] = draw;
        }
        full_grid.GlobalSumVector(draws.data(), replica_count);
        assert(draws[0] != draws[1]);

        LatticeGaugeField U(&split_grid);
        SU<Nc>::ColdConfiguration(parallel_rng, U);

        WilsonGaugeActionR bare_action(1.0);
        ConstrainedAction<WilsonGaugeActionR> action(bare_action, {1.0, 10.0 * replica_index, 100.0});
        typedef RobbinsMonroSolver<ConstrainedAction<WilsonGaugeActionR> > Solver;
        Solver solver(action, RobbinsMonroParameters(0, 0, 10, 1.0));
        LLRReplicaExchangeModule<Solver> module(full_grid, split_grid, replica_index, replica_index, solver, 1, exchange_rng);

        module.TrajectoryComplete(1, U, serial_rng, parallel_rng);
        assert(module.exchange_attempt_count() == 1);

        // Check 2: check that the intervals have correctly swapped.
        int expected_interval = 1 - replica_index;
        assert(module.interval_index() == expected_interval);
        uint64_t interval_index_sum = module.interval_index();
        split_grid.GlobalSum(interval_index_sum);
        assert(interval_index_sum == static_cast<uint64_t>(expected_interval) * split_grid.ProcessorCount());

        // Check 3: Check that the action was correctly updated during the interval swap.
        assert(action.parameters().S0 == 10.0 * expected_interval);
        assert(solver.status().trajectories_remaining_in_phase == 9);

        // Check 4: Check that the accept/reject state is correctly
        // known to the full grid.
        RealD exchange_draw = 0.0;
        {
            ScopedWorldCommunicator exchange_scope(full_grid);
            random(exchange_rng, exchange_draw);
        }
        RealD exchange_draw_sum = exchange_draw;
        full_grid.GlobalSum(exchange_draw_sum);
        assert(std::abs(exchange_draw_sum - exchange_draw * full_grid.ProcessorCount()) < 1.0e-12);
    }

    Grid_finalize();
    return 0;
}

