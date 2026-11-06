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
#pragma once

#include <Grid/Grid.h>
#include <Grid/qcd/llr/LLRSolverModule.h>

using namespace Grid;

typedef ConstrainedAction<WilsonGaugeActionR> ConstrainedWilsonGaugeAction;

template<class Solver>
bool run_test(ConstrainedActionParameters action_parameters,
	      LLRSolverParameters solver_parameters,
	      std::string test_name, int n_traj, RealD a_expected, RealD epsilon_a)
{
  RealD S0 = action_parameters.S0;
  RealD a0 = action_parameters.a;
  RealD sigma = action_parameters.sigma;

  // The HMC runner
  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;
  HMCWrapper TheHMC;

  // The grid
  TheHMC.Resources.AddFourDimGrid("gauge");
  
  // Checkpointer definition
  CheckpointerParameters CPparams;  
  CPparams.config_prefix = "ckpoint_lat";
  CPparams.rng_prefix = "ckpoint_rng";
  CPparams.saveInterval = 65535;
  CPparams.format = "IEEE64BIG";
  
  TheHMC.Resources.LoadNerscCheckpointer(CPparams);
  
  // The RNG
  RNGModuleParameters RNGpar;
  RNGpar.serial_seeds = "1 2 3 4 5";
  RNGpar.parallel_seeds = "6 7 8 9 10";
  TheHMC.Resources.SetRNGSeeds(RNGpar);
  
  // Construct observables
  typedef PlaquetteMod<HMCWrapper::ImplPolicy> PlaqObs;
  TheHMC.Resources.AddObservable<PlaqObs>();
  
  typedef TopologicalChargeMod<HMCWrapper::ImplPolicy> QObs;
  TopologyObsParameters TopParams;
  TopParams.interval = 1;
  TopParams.do_smearing = false;
  TopParams.Smearing.init_step_size = 0.01;
  TopParams.Smearing.tolerance = 1e-5;
  //  TopParams.Smearing.steps = 200;
  //  TopParams.Smearing.step_size = 0.01;
  TopParams.Smearing.meas_interval = 50;
  TopParams.Smearing.maxTau = 2.0; 
  TheHMC.Resources.AddObservable<QObs>(TopParams);
  
  // LLR updater
  WilsonGaugeActionR bare_action(1.0); // beta HAS to be 1 for LLR
  ConstrainedWilsonGaugeAction constrained_action(bare_action, action_parameters);
  Solver solver(constrained_action, solver_parameters);
  typedef LLRSolverModule<Solver> LLRMod;
  TheHMC.Resources.AddObservable<LLRMod>(solver);
  TheHMC.Resources.AddObservableName("LLRSolver");
        
  // Collect actions
  ActionLevel<HMCWrapper::Field> Level1(1);
  Level1.push_back(&constrained_action);
  TheHMC.TheAction.push_back(Level1);
  
  // HMC parameters are serialisable
  TheHMC.Parameters.Trajectories = n_traj;
  TheHMC.Parameters.NoMetropolisUntil = 0;
  TheHMC.Parameters.MD.MDsteps = 100;
  TheHMC.Parameters.MD.trajL   = 1.0;

  TheHMC.Run();  // no smearing

  // Report intermediate action values
  auto theObservables = TheHMC.Resources.GetObservables();
  int iObs = TheHMC.Resources.GetObservableIndex("LLRSolver");
  typedef LLRActionLogger<Solver> ActionLogger;
  ActionLogger* theLLRmodule = reinterpret_cast<ActionLogger*>(theObservables[iObs]);
  std::vector<std::pair<int, RealD> > theActions = theLLRmodule->getMeanActionVector();
  std::vector<std::pair<int, RealD> > theAs = theLLRmodule->getAVector();
  std::cout << GridLogMessage << "Test_" << test_name << "_solver: Pairs of trajectory-mean_action for every iteration: " << theActions << std::endl;
    std::cout << GridLogMessage << "Test_" << test_name << "_solver: Pairs of trajectory-a for every iteration: " << theAs << std::endl;

  // test final values
  RealD mean_action = solver.status().last_update.mean_action;
  bool S0_converged = std::abs(mean_action - S0) < 2 * sigma;
  bool a_converged  = std::abs(constrained_action.parameters().a - a_expected) < epsilon_a;
  std::cout << GridLogMessage << "Test_" << test_name << "_solver: " << " a_initial = " << a0 <<
    ", S0 = " << S0 << (S0_converged ? ", S0 converged" : ", S0 did not converge") <<
    (a_converged ? ", a converged" : ", a did not converge") << std::endl;
  std::cout << GridLogMessage << "Test_" << test_name << "_solver: a_final = " << constrained_action.parameters().a <<
    ", Sunconstrained = " << mean_action << std::endl;

  return S0_converged && a_converged;
}

