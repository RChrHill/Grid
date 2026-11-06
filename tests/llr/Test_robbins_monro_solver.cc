/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./tests/Test_llr_robbinsmonro.cc

Copyright (C) 2015, 2026

Author: Peter Boyle <pabobyle@ph.ed.ac.uk>
Author: neo <cossu@post.kek.jp>
Author: Guido Cossu <guido.cossu@ed.ac.uk>
Author: Ryan Hill <guido.cossu@ed.ac.uk>
Author: Ilektra Christidi <ilektra.christidi@ucl.ac.uk>

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
#include <Grid/Grid.h>
#include <Grid/qcd/llr/LLRSolverModule.h>
#include "run_llr_test.cc"

// TODO: unit tests for edge cases of number of trajectories in various phases (eg. rethrm=0, therm=0, correct phase at end?)

int main(int argc, char **argv) 
{
  using namespace Grid;

  // Ignore command line arguments and enforce --grid 4.4.4.4
  int updatedArgc = 3;
  char** updatedArgv = (char**)malloc(sizeof(char*) * 4);
  updatedArgv[0] = argv[0];
  updatedArgv[1] = (char*)(new std::string("--grid"))->c_str();
  updatedArgv[2] = (char*)(new std::string("4.4.4.4"))->c_str();
  updatedArgv[3] = nullptr;
  Grid_init(&updatedArgc, &updatedArgv);
  GridLogLayout();

  // Starting and expected final values for constrained action and solver parameters
  RealD a0, S0, a_expected, sigma  = 3.0;
  int n_thermalisation = 100, n_rethermalisation = 10, n_accumulation = 5, n_iter = 14;
  int n_traj = n_thermalisation + n_accumulation + (n_iter -1) * (n_accumulation + n_rethermalisation); // 300 for default values
  RealD epsilon_a = 0.3;
  int n_steps = 2;
  std::vector<RealD> a0_scan(n_steps);
  std::vector<RealD> a_expect_scan(n_steps);
  std::vector<RealD> S0_scan(n_steps);
  S0_scan = {768, 870};
  std::string gaugeGroup = "SU(3)";
  if (Sp2n_config) {
    a0_scan = {-10, 30};
    a_expect_scan = {11, 12};
    gaugeGroup = "Sp(4)";
  }
  else {
    a0_scan = {-10, 20};
    a_expect_scan = {5.5, 5.7};
  }

  // Scan starting values for a (including negative and extremes) and S0 and run test
  std::vector<bool> converged(n_steps*n_steps);
  for (int i=0; i<n_steps; i++) {
    S0 = S0_scan[i];
    a_expected = a_expect_scan[i];
    std::cout << GridLogMessage << "Test_robbins_monro_solver: Comparing with expected values for " <<
      gaugeGroup << ": S0 = " << S0 << ", a = " << a_expected << std::endl;
    
    for (int j=0; j<n_steps; j++) {
      a0 = a0_scan[j];

      typedef LLRSolver<ConstrainedWilsonGaugeAction, &RobbinsMonro> Solver;
      ConstrainedActionParameters action_parameters{.a=a0, .S0=S0, .sigma=sigma};
      LLRSolverParameters solver_parameters(n_thermalisation, n_rethermalisation, n_accumulation, 1.0);
      converged[i*n_steps+j] = run_test<Solver>(action_parameters, solver_parameters, "robbins_monro", n_traj, a_expected, epsilon_a);
    }
  }
  
  assert(std::all_of(converged.begin(), converged.end(), [](bool v) { return v; }));

  Grid_finalize();

} // main
