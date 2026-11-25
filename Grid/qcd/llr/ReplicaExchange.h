/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ReplicaExchange.h

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

#include <Grid/qcd/action/gauge/ConstrainedAction.h>

NAMESPACE_BEGIN(Grid);

// Result of an attempted interval swap (first_interval_index, first_interval_index + 1)
struct ReplicaExchangeDecision
{
    int first_interval_index;
    int first_replica_index;
    int second_replica_index;
    RealD log_acceptance;
    bool accepted;
};

// Info given by replicas to facilitate interval exchange
struct ReplicaSample
{
    int interval_index;
    ConstrainedActionParameters action_parameters;
    RealD bare_action;
};

// ln P = U_1(S_1) + U_2(S_2) - U_1(S_2) - U_2(S_1) for swapping the intervals held by two replicas
inline RealD logExchangeAcceptance(const ReplicaSample &first, const ReplicaSample &second)
{
    return first.action_parameters.constrain(first.bare_action)
         + second.action_parameters.constrain(second.bare_action)
         - first.action_parameters.constrain(second.bare_action)
         - second.action_parameters.constrain(first.bare_action);
}

// Loop over either odd or even intervals, and determine which are to be exchanged
// based on an accept/reject.
inline std::vector<ReplicaExchangeDecision> decideAcceptRejects(unsigned int exchange_attempt, const std::vector<ReplicaSample> &samples, const std::vector<RealD>& metropolis_draws)
{
    const int replica_count = samples.size();
    assert(metropolis_draws.size() >= replica_count / 2);

    // Generate a map of interval index -> replica index.
    // This could be stored on a class as a member, with `decideAcceptRejects` as a member fn.
    std::vector<int> interval_to_replica(replica_count, -1);
    for (int replica = 0; replica < replica_count; ++replica)
    {
        assert(interval_to_replica.at(samples[replica].interval_index) == -1);
        interval_to_replica[samples[replica].interval_index] = replica;
    }

    // For each interval covered by this attempt, take the RNG draw for that interval and convert that
    // into an accept/reject bool.
    // Also log which replicas etc. that swap attempt corresponds to.
    std::vector<ReplicaExchangeDecision> decisions;
    for (int interval = exchange_attempt % 2; interval + 1 < replica_count; interval += 2)
    {
        ReplicaExchangeDecision decision;
        decision.first_interval_index = interval;
        decision.first_replica_index  = interval_to_replica[interval];
        decision.second_replica_index = interval_to_replica[interval + 1];
        assert(samples[decision.first_replica_index].action_parameters.S0 < samples[decision.second_replica_index].action_parameters.S0);
        decision.log_acceptance       = logExchangeAcceptance(samples[decision.first_replica_index], samples[decision.second_replica_index]);
        decision.accepted             = std::log(metropolis_draws[interval / 2]) < decision.log_acceptance;
        decisions.push_back(decision);
    }
    return decisions;
}

NAMESPACE_END(Grid);
