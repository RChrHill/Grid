/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ReplicaExchangeModule.h

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

#include <Grid/qcd/llr/ReplicaExchange.h>
#include <Grid/qcd/llr/RobbinsMonroSolverModule.h>
#include <Grid/qcd/llr/ScopedWorldCommunicator.h>

NAMESPACE_BEGIN(Grid);

struct IntervalState
{
    int interval_index;
    ConstrainedActionParameters action;
    RobbinsMonroStatus solver;
};

template <class SolverType>
class LLRReplicaExchangeModule
    : public HMCModuleBase<HmcObservable<typename SolverType::Field> >,
      public HmcObservable<typename SolverType::Field>
{
public:
    typedef typename SolverType::Field Field;
    typedef HMCModuleBase<HmcObservable<Field> > Base;
    typedef typename Base::Product Product;
private:
    // Everything a replica contributes to an exchange step.
    struct Record
    {
        ReplicaSample sample;
        IntervalState interval;
        int previous_iteration;
    };
public:
    LLRReplicaExchangeModule(GridBase &full_grid, GridBase &split_grid, int replica_index, int interval_index,
                             SolverType &solver, int steps_between_exchanges, GridSerialRNG &exchange_rng)
        : exchange_rng_(exchange_rng), solver_(solver), full_grid_(full_grid), split_grid_(split_grid),
          replica_index_(replica_index), interval_index_(interval_index), steps_between_exchanges_(steps_between_exchanges),
          replica_count_(full_grid.ProcessorCount() / split_grid.ProcessorCount()), exchange_attempt_count_(0)
    {
        assert(steps_between_exchanges_ > 0);
        assert(replica_index_ >= 0 && replica_index_ < replica_count_);
        assert(interval_index_ >= 0 && interval_index_ < replica_count_);

        
        // Copy the solver parameters to all ranks and sanity-check that they were copied.
        const RobbinsMonroParameters& local = solver_.parameters();
        for (const RobbinsMonroParameters& other : duplicateFromBossToAllRanksInReplica(local))
        {
            assert(other.initial_thermalisation_trajectories == local.initial_thermalisation_trajectories);
            assert(other.rethermalisation_trajectories == local.rethermalisation_trajectories);
            assert(other.trajectories_per_update == local.trajectories_per_update);
            assert(other.gain == local.gain);
        }
    }

    virtual Product *getPtr()
    {
        return this;
    }

    virtual void TrajectoryComplete(int trajectory, Field &U, GridSerialRNG &, GridParallelRNG &)
    {
        complete_trajectory(trajectory, U);
    }

    virtual void TrajectoryComplete(int trajectory, ConfigurationBase<Field> &configuration, GridSerialRNG &, GridParallelRNG &, bool)
    {
        Field &U = configuration.get_U(solver_.action_.is_smeared);
        complete_trajectory(trajectory, U);
    }

    int interval_index() const
    {
        return interval_index_;
    }

    unsigned int exchange_attempt_count() const
    {
        return exchange_attempt_count_;
    }
private:
    void complete_trajectory(int trajectory, Field &U)
    {
        int previous_iteration = solver_.status().iteration;
        if (trajectory % steps_between_exchanges_ != 0)
        {
            solver_.record_configuration(trajectory, U);
            report_local_update(previous_iteration);
            return;
        }

        Record local;
        local.sample = {interval_index_,solver_.action_.parameters(), solver_.action_.Sunconstrained(U)};
        solver_.record_configuration(trajectory, U);
        local.interval = interval_state();
        local.previous_iteration = previous_iteration;

        std::vector<Record> records = duplicateFromBossToAllRanksInReplica(local);
        std::vector<ReplicaSample> samples;
        for (const Record &record : records)
        {
            samples.push_back(record.sample);
        }

        // GridSerialRNG::fill broadcasts over the world communicator, here we're
        // swapping out the replica communicator for the world communicator so the
        // same value can be distributed across the entire grid.
        std::vector<RealD> metropolis_draws(replica_count_ / 2);
        {
            ScopedWorldCommunicator scope(full_grid_);
            for (RealD &draw : metropolis_draws)
            {
                random(exchange_rng_, draw);
            }
        }

        // Convert the random numbers into accept/reject for each potential interval exchange.
        std::vector<ReplicaExchangeDecision> decisions = decideAcceptRejects(exchange_attempt_count_, samples, metropolis_draws);
        if (full_grid_.IsBoss())
        {
            report_exchanges(decisions, trajectory, exchange_attempt_count_);
            report_updates(records);
        }

        // Perform the accepted interval exchanges.
        for (const ReplicaExchangeDecision &decision : decisions)
        {
            if (decision.accepted)
            {
                std::swap(records[decision.first_replica_index].interval, records[decision.second_replica_index].interval);
            }
        }
        restore_interval_state(records[replica_index_].interval);
        ++exchange_attempt_count_;
    }

    IntervalState interval_state() const
    {
        IntervalState state;
        state.interval_index = interval_index_;
        state.action = solver_.action_.parameters();
        state.solver = solver_.status();
        return state;
    }

    void restore_interval_state(const IntervalState &state)
    {
        assert(state.interval_index >= 0 && state.interval_index < replica_count_);
        solver_.action_.set_parameters(state.action);
        solver_.restore_state(state.solver);
        interval_index_ = state.interval_index;
    }

    // Every rank in a replica gets the value on the boss rank
    template <class T>
    std::vector<T> duplicateFromBossToAllRanksInReplica(const T& local)
    {
        std::vector<T> gathered(replica_count_);
        std::memset(gathered.data(), 0, gathered.size() * sizeof(T));   // <¬
        if (split_grid_.IsBoss())                                       //  | 
        {                                                               //  |
            std::memcpy(&gathered[replica_index_], &local, sizeof(T));  //  |
        }   //          ^                                               //  |
            //          |                                               //  |
            //       /-------------------\                                  |
        // All ranks except the boss value are set to 0                     |
        //                                     \------/---------------------/
        //
        // GlobalSumVector therefore adds the boss rank value to a bunch of 0s, and copies
        // that to all ranks, so all ranks get the boss value.
        full_grid_.GlobalSumVector(reinterpret_cast<uint64_t *>(gathered.data()), static_cast<int>(gathered.size() * sizeof(T) / sizeof(uint64_t)));
        return gathered;
    }

    void report_local_update(int previous_iteration)
    {
        RobbinsMonroStatus status = solver_.status();
        if (split_grid_.IsBoss() && status.iteration != previous_iteration)
        {
            LogRobbinsMonroUpdate(status.last_update, " interval " + std::to_string(interval_index_));
        }
    }

    static void report_updates(const std::vector<Record> &records)
    {
        for (const Record &record : records)
        {
            if (record.interval.solver.iteration != record.previous_iteration)
            {
                LogRobbinsMonroUpdate(record.interval.solver.last_update, " interval " + std::to_string(record.interval.interval_index));
            }
        }
    }

    static void report_exchanges(const std::vector<ReplicaExchangeDecision> &decisions, unsigned int trajectory, unsigned int exchange_attempt)
    {
        for (const ReplicaExchangeDecision &decision : decisions)
        {
            std::cout << GridLogMessage        << "[ReplicaExchange] trajectory " << trajectory
                      << " attempt "           << exchange_attempt
                      << " intervals "         << decision.first_interval_index
                      << " and "               << decision.first_interval_index + 1
                      << " replicas "          << decision.first_replica_index
                      << " and "               << decision.second_replica_index
                      << " log acceptance "    << decision.log_acceptance
                      << " exchange accepted " << decision.accepted << std::endl;
        }
    }
private:
    GridSerialRNG &exchange_rng_;
    SolverType &solver_;
    GridBase &full_grid_;
    GridBase &split_grid_;
    int replica_index_;
    int interval_index_;
    int steps_between_exchanges_;
    int replica_count_;
    unsigned int exchange_attempt_count_;
};

NAMESPACE_END(Grid);
