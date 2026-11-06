/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: ./examples/Example_LLR.cc

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

using namespace Grid;

/*! @brief Example HMC run of the LLR workflow.
 *
 * Uses the constrained action
 *   \f$\beta S[U] \rightarrow aS[U] + (S[U]-S_0)^2/(2\sigma^2)\f$
 *
 * Parameters that can be provided in the command line:
 * - the initial \f$a\f$: `--a_init <value>`. Default = 1
 * - the initial \f$S_0\f$: `--S0_init <value>`. Default = 1000
 * - the LLR solver parameters and number of iterations:
 *   `--llr_params <n_thermalisation>,<n_rethermalisation>,<n_accumulation>,<n_iter>`.
 *   Default = (100, 10, 5, 14)
 *   The number of trajectories is calculated from those parameters.
 *
 * The LLR algorithm can be turned off with `--noLLR`, in which case
 * the unconstrained action is used.
 * Parameters that can be provided in the command line:
 * - the \f$beta\f$ of the unconstrained action: `--beta <value>`. Default = 1
 * - the number of trajectories: `--n_traj <value>`. Default = 300
 */


// User-defined updater - it can be anything, this is just an example
void NullUpdater(LLRSolverUpdate &update, RealD gain, RealD sigma, int iteration)
{
  std::cout << GridLogMessage << "Example_LLR: Calling custom updater." << std::endl;
  update.updated_a = update.previous_a;
}


int main(int argc, char **argv) 
{
  using namespace Grid;

  Grid_init(&argc, &argv);
  GridLogLayout();

  std::string gaugeGroup = "SU(3)";
  if (Sp2n_config) {
    gaugeGroup = "Sp(4)";
  }

  // Parse command line for action and solver parameters
  bool doLLR = true;
  if( GridCmdOptionExists(argv, argv+argc, "--noLLR") ) {
    std::cout << GridLogMessage << "Example_LLR: no-LLR for " << gaugeGroup << ": turning off the constrained action and LLR solver. Arguments --a_init, --S0_init, and --llr_params will be ignored, but --beta and --n_traj can be provided." << std::endl;
    doLLR = false;
  }
  else {
    std::cout << GridLogMessage << "Example_LLR: LLR for " << gaugeGroup << ": running the constrained action and LLR solver. Argument --beta and --n_traj will be ignored, but --a_init, --S0_init, and --llr_params can be provided." << std::endl;
  }
  RealD alpha0 = 1, S0  = 1000, sigma  = 3.0;
  RealD beta = 1; // HAS to be 1 for LLR
  int n_thermalisation = 100, n_rethermalisation = 10, n_accumulation = 5, n_iter = 14;
  int n_traj = 300;
  std::string arg;
  std::vector<int> vec;
  if (doLLR) {
    if( GridCmdOptionExists(argv, argv+argc, "--a_init") ) {
      arg = GridCmdOptionPayload(argv, argv+argc, "--a_init");
      GridCmdOptionFloat(arg, alpha0);
    }
    if( GridCmdOptionExists(argv, argv+argc, "--S0_init") ) {
      arg = GridCmdOptionPayload(argv, argv+argc, "--S0_init");
      GridCmdOptionFloat(arg, S0);
    }
    if( GridCmdOptionExists(argv, argv+argc, "--llr_params") ) {
       arg = GridCmdOptionPayload(argv, argv+argc, "--llr_params");
       GridCmdOptionIntVector(arg,vec);
       n_thermalisation = vec[0];
       n_rethermalisation = vec[1];
       n_accumulation = vec[2];
       n_iter = vec[3];
       n_traj = n_thermalisation + n_accumulation + (n_iter -1) * (n_accumulation + n_rethermalisation);
    }
    std::cout << GridLogMessage << "Example_LLR: constrained action parameters: S0 = " << S0 << ", a = " << alpha0 << ", sigma = " << sigma << std::endl;
    std::cout << GridLogMessage << "Example_LLR: LLR solver parameters: n_thermalisation = " << n_thermalisation << ", n_rethermalisation = " << n_rethermalisation << ", n_accumulation = " << n_accumulation << ", n_iter = " << n_iter << ", n_traj = " << n_traj << std::endl;
  }
  else {
    if( GridCmdOptionExists(argv, argv+argc, "--beta") ) {
      arg = GridCmdOptionPayload(argv, argv+argc, "--beta");
      GridCmdOptionFloat(arg, beta);
    }
    if( GridCmdOptionExists(argv, argv+argc, "--n_traj") ) {
      arg = GridCmdOptionPayload(argv, argv+argc, "--n_traj");
      GridCmdOptionInt(arg, n_traj);
    }
    std::cout << GridLogMessage << "Example_LLR: unconstrained action beta = " << beta << std::endl;
    std::cout << GridLogMessage << "Example_LLR: MD trajectories n_traj = " << n_traj << std::endl;
  }
  
  // The HMC runner
  typedef GenericHMCRunner<MinimumNorm2> HMCWrapper;  // Uses the default minimum norm
  HMCWrapper TheHMC;
  
  // Grid from the command line
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
  
  WilsonGaugeActionR bare_action(beta);
  typedef ConstrainedAction<WilsonGaugeActionR> ConstrainedWilsonGaugeAction;
  ConstrainedActionParameters action_parameters{.a=alpha0, .S0=S0, .sigma=sigma};
  ConstrainedWilsonGaugeAction constrained_action(bare_action, action_parameters);
  // Null updater (user-defined)
  typedef LLRSolver<ConstrainedWilsonGaugeAction, &NullUpdater> SolverN;
  SolverN solverN(constrained_action, LLRSolverParameters(n_thermalisation, n_rethermalisation, n_accumulation, 1.0));
  typedef LLRSolverModule<SolverN> NullMod;
  // Robbins-Monro updater
  typedef LLRSolver<ConstrainedWilsonGaugeAction, &RobbinsMonro> Solver;
  Solver solver(constrained_action, LLRSolverParameters(n_thermalisation, n_rethermalisation, n_accumulation, 1.0));
  typedef LLRSolverModule<Solver> RmMod;
  if (doLLR)
    {
      TheHMC.Resources.AddObservable<NullMod>(solverN);
      TheHMC.Resources.AddObservable<RmMod>(solver);
    }
  
  // Collect actions
  ActionLevel<HMCWrapper::Field> Level1(1);
  if (doLLR)
    Level1.push_back(&constrained_action);
  else
    Level1.push_back(&bare_action);    
  TheHMC.TheAction.push_back(Level1);
    
  // HMC parameters are serialisable
  TheHMC.Parameters.Trajectories = n_traj;
  TheHMC.Parameters.NoMetropolisUntil = 0;
  TheHMC.Parameters.MD.MDsteps = 100;
  TheHMC.Parameters.MD.trajL   = 1.0;
  
  TheHMC.ReadCommandLine(argc, argv); // these can be parameters from file
  
  TheHMC.Run();  // no smearing

  // Report final values
  if (doLLR) {
    RealD mean_actionN = solverN.status().last_update.mean_action;
    RealD mean_action = solver.status().last_update.mean_action;
    RealD last_aN = solverN.status().last_update.updated_a;
    RealD last_a = solver.status().last_update.updated_a;
    RealD a_final  = constrained_action.parameters().a;
    std::cout << GridLogMessage << "Example_LLR: a_final = " << a_final << ", last_a = " << last_a << ", last_a Null = " << last_aN << ", Sunconstrained = " << mean_action << ", Sunconstrained Null = " << mean_actionN << std::endl;
  }
  
  Grid_finalize();

} // main
