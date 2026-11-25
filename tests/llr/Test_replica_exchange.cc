/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Test_replica_exchange.cc

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


// Just a plain test of the replica exchange correctness against a reference probability

#include <Grid/Grid.h>
#include <Grid/qcd/llr/ReplicaExchange.h>

using namespace Grid;

namespace
{

bool approximately_equal(RealD actual, RealD expected)
{
    return std::abs(actual - expected) <= 1.0e-12 * (1.0 + std::abs(expected));
}

// A sample's replica is its position in the vector passed to decideAcceptRejects
ReplicaSample sample(int interval_index, RealD a, RealD S0, RealD sigma, RealD bare_action)
{
    ReplicaSample result;
    result.interval_index    = interval_index;
    result.action_parameters = {a, S0, sigma};
    result.bare_action       = bare_action;
    return result;
}

// Checks acceptance with hardcoded probability inputs
void test_acceptance()
{
    ReplicaSample first  = sample(0, 1.0, 0.0, 2.0, 1.0);
    ReplicaSample second = sample(1, 1.5, 4.0, 2.0, 5.0);

    assert(approximately_equal(logExchangeAcceptance(first, second), -2.0));
    assert(approximately_equal(logExchangeAcceptance(second, first), -2.0));

    std::vector<ReplicaSample> pair = {first, second};
    std::vector<ReplicaExchangeDecision> accepted = decideAcceptRejects(0, pair, {0.1});
    assert(accepted.size() == 1);
    assert(accepted[0].accepted);
    assert(accepted[0].first_replica_index  == 0);
    assert(accepted[0].second_replica_index == 1);

    std::vector<ReplicaExchangeDecision> rejected = decideAcceptRejects(0, pair, {0.2});
    assert(!rejected[0].accepted);
}

// Checks that exchanges are based on interval neighbours rather than replica neighbours
void test_interval_pairing()
{
    std::vector<ReplicaSample> samples;
    samples.push_back(sample(2, 0.0, 2.0, 1.0, 0.0)); // replica 0
    samples.push_back(sample(0, 0.0, 0.0, 1.0, 0.0)); // replica 1
    samples.push_back(sample(3, 0.0, 3.0, 1.0, 0.0)); // replica 2
    samples.push_back(sample(1, 0.0, 1.0, 1.0, 0.0)); // replica 3

    std::vector<ReplicaExchangeDecision> even = decideAcceptRejects(0, samples, {0.5, 0.5});
    assert(even.size() == 2);
    assert(even[0].first_replica_index == 1 && even[0].second_replica_index == 3);
    assert(even[1].first_replica_index == 0 && even[1].second_replica_index == 2);

    std::vector<ReplicaExchangeDecision> odd = decideAcceptRejects(1, samples, {0.5, 0.5});
    assert(odd.size() == 1);
    assert(odd[0].first_interval_index == 1);
    assert(odd[0].first_replica_index == 3 && odd[0].second_replica_index == 0);
}

}  // namespace

int main(int argc, char **argv)
{
    Grid_init(&argc, &argv);

    test_acceptance();
    test_interval_pairing();

    Grid_finalize();
    return 0;
}

