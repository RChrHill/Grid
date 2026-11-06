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
  RealD a0, S0 = 768, a_expected, sigma  = 3.0;
  int n_thermalisation = 250, n_rethermalisation = 20, n_accumulation = 5, n_iter = 190;
  int n_traj = n_thermalisation + n_accumulation + (n_iter -1) * (n_accumulation + n_rethermalisation);
  RealD epsilon_a = 0.03;
  std::string gaugeGroup = "SU(3)";
  if (Sp2n_config) {
    a0 = 9;
    a_expected = 9;
    gaugeGroup = "Sp(4)";
  }
  else {
    a0 = 5.5;
    a_expected = 5.5;
  }

  // Run test
  bool converged;
  std::cout << GridLogMessage << "Test_robbins_monro_solver: Comparing with expected values for " <<
    gaugeGroup << ": S0 = " << S0 << ", a = " << a_expected << std::endl;
    
  typedef LLRSolver<ConstrainedWilsonGaugeAction, &RobbinsMonro> Solver;
  ConstrainedActionParameters action_parameters{.a=a0, .S0=S0, .sigma=sigma};
  LLRSolverParameters solver_parameters(n_thermalisation, n_rethermalisation, n_accumulation, 1.0);
  converged = run_test<Solver>(action_parameters, solver_parameters, "robbins_monro", n_traj, a_expected, epsilon_a);
  
  assert(converged);

  Grid_finalize();

} // main
