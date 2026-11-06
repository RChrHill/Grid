/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: LLRSolver.h

Copyright (C) 2026

Author: Ryan Hill <Ryan.Hill@ed.ac.uk>
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

#include <cassert>

#include <Grid/qcd/action/gauge/ConstrainedAction.h>

/*! \file
 *
 */
/// \cond DO_NOT_DOCUMENT
NAMESPACE_BEGIN(Grid);
/// \endcond

enum class LLRSolverPhase
/*! @brief The LLR solver phase.
 *
 * It can either be Thermalising or Accumulating. When/How long the solver
 * stays in each phase is controlled by the LLRSolverParameters structure.
 */
{
  Thermalising,
  Accumulating
};

struct LLRSolverParameters
/*! @brief The parameters controlling the LLR solver.
 *
 * Those include, among others, the number of trajectories in each phase.
 */
{
  int initial_thermalisation_trajectories; ///< @brief The number of trajectories to use for initial thermalisation
  int rethermalisation_trajectories; ///< @brief The number of trajectories to use for rethermalisation, after the accumulation phase
  int trajectories_per_update; ///< @brief The number of trajectories to accumulate per update
  RealD gain; ///< @brief The gain used to calculate the new \f$a\f$ in every update

  /*! @brief Construct a LLRSolverParameters structure.
   * @param[in] initial_thermalisation_trajectories_: the number of trajectories
   *            to use for initial thermalisation
   * @param[in] rethermalisation_trajectories_: the number of trajectories
   *            to use for rethermalisation, after the accumulation phase
   * @param[in] trajectories_per_update_: the number of trajectories to accumulate
   *            per update
   * @param[in] gain_: the gain used to calculate the new \f$a\f$ in every update
   */
  LLRSolverParameters(int initial_thermalisation_trajectories_ = 0,
                         int rethermalisation_trajectories_ = 0,
                         int trajectories_per_update_ = 1,
                         RealD gain_ = 1.0)
      : initial_thermalisation_trajectories(initial_thermalisation_trajectories_),
        rethermalisation_trajectories(rethermalisation_trajectories_),
        trajectories_per_update(trajectories_per_update_),
        gain(gain_)
  {}
};

struct LLRSolverUpdate
/*! @brief Structure holding useful information every time \f$a\f$ is updated
 * by the LLR solver.
 */
{
  int trajectory; ///< @brief The current trajectory at update time
  int iteration; ///< @brief The current iteration at update time
  RealD mean_action; ///< @brief The mean unconstrained action per trajectory in this update
  RealD residual; ///< @brief The difference between mean_action and S0
  RealD previous_a; ///< @brief The \f$a\f$ parameter of the previous update
  RealD updated_a; ///< @brief The \f$a\f$ parameter of the current update
};

struct LLRSolverStatus
/*! @brief Structure holding the LLR solver status.
 *
 * Needed for replica-swapping.
 */
{
  LLRSolverPhase phase; ///< @brief The current phase of the solver (Thermalising or Accumulating)
  int iteration; ///< @brief The current iteration
  int trajectories_remaining_in_phase; ///< @brief The number of trajectories remaining in the current phase
  RealD accumulated_action; ///< @brief The unconstrained action accumulated so far in this update interval
  LLRSolverUpdate last_update; ///< @brief The last \f$a\f$ update information structure
};

/*! @brief Predefined solver updaters (currently Robbins-Monro and Newton-Raphson).
 * The user can also define their own locally in their application.
 *
 * @param[inout] update: The LLRSolverUpdate structure to modify
 * @param[in] gain: the gain used to calculate the new \f$a\f$ in every update
 * @param[in] sigma: the constrained action sigma
 * @param[in] iteration: the current iteration at update time
 */
typedef void (*updater)(LLRSolverUpdate&, RealD, RealD, int);
void RobbinsMonro(LLRSolverUpdate &update, RealD gain, RealD sigma, int iteration)
{
  update.updated_a = update.previous_a + gain * update.residual / (sigma * sigma * iteration);
}

void NewtonRaphson(LLRSolverUpdate &update, RealD gain, RealD sigma, int iteration)
{
  update.updated_a = update.previous_a + gain * update.residual / (sigma * sigma);
}


template <class ConstrainedActionType, updater Updater>
class LLRSolver
/*! @brief The LLR solver.
 *
 * Expects two template parameters:
 * @param ConstrainedActionType: The constrained action type.
 * @param Updater: The solver updater.
 */
{
public:
  typedef ConstrainedActionType ActionType;
  typedef typename ActionType::GaugeField Field;
  
  /*! @brief Construct an LLR solver.
   * @param[in] action: the constrained action
   * @param[in] parameters: the LLRSolverParameters parameter structure
   */
  LLRSolver(ActionType &action, LLRSolverParameters parameters)
      : action_(action)
      , parameters_(parameters)
      , status_{
          .phase=LLRSolverPhase::Accumulating,
          .iteration=1,
          .trajectories_remaining_in_phase=parameters.trajectories_per_update,
          .accumulated_action=0
        }
  {
    if (parameters_.initial_thermalisation_trajectories != 0)
    {
      status_.phase = LLRSolverPhase::Thermalising;
      status_.trajectories_remaining_in_phase = parameters_.initial_thermalisation_trajectories;
    }
  }

  /*! @brief Call the appropriate phase of the LLR solver based on
   * the current status.
   * @param[in] trajectory: the current trajectory number
   * @param[in] U: the gauge field
   */
  void record_configuration(int trajectory, Field &U)
  {
    status_.trajectories_remaining_in_phase--;
    if (status_.phase == LLRSolverPhase::Accumulating)
    {
      std::cout << GridLogMessage << "LLR Phase: Accumulating (remaining trajectories: " << status_.trajectories_remaining_in_phase << ")" << std::endl;
      accumulate(trajectory, U);
    }
    else // Thermalising
    {
      std::cout << GridLogMessage << "LLR Phase: Thermalising (remaining trajectories: " << status_.trajectories_remaining_in_phase << ")" << std::endl;
      thermalise();
    }
  }

  /*! @brief Call record_configuration for the current trajectory and field configuration.
   * @param[in] trajectory: the current trajectory number
   * @param[in] configuration: the gauge field configuration
   */
  void record_configuration(int trajectory, ConfigurationBase<Field> &configuration)
  {
    Field &U = configuration.get_U(action_.is_smeared);
    record_configuration(trajectory, U);
  }

  // Accessors for replica-swapping
  /*! @brief Getter for the LLR solver parameters.
   * @returns The solver parameters LLRSolverParameters structure
   */
  const LLRSolverParameters &parameters() const { return parameters_; }
  /*! @brief Restore the LLR solver state.
   * @param[in] state: the LLRSolverStatus status structure
   */
  void restore_state(const LLRSolverStatus &state) { status_ = state; }
  /*! @brief Getter for the LLR solver status object.
   * @returns The solver LLRSolverStatus status structure
   */  
  LLRSolverStatus status() const { return status_; }
  
public:
  ActionType &action_; ///< @brief The constrained action
private:
  LLRSolverParameters parameters_; ///< @brief The solver parameters structure
  LLRSolverStatus status_; ///< @brief The solver status structure

  /*! @brief Call the accumulate phase of the LLR solver for
   * the current trajectory.
   * @param[in] trajectory: the current trajectory number
   * @param[in] U: the gauge field
   */
  void accumulate(int trajectory, Field &U)
  {
    status_.accumulated_action += action_.Sunconstrained(U);
    if (status_.trajectories_remaining_in_phase > 0) // Ready to update? If not, return to continue sampling.
    {
      return;
    }

    // After we've gathered enough samples, perform the update of a.
    LLRSolverUpdate update;
    update.trajectory = trajectory;
    update.iteration = status_.iteration;
    update.mean_action = status_.accumulated_action / parameters_.trajectories_per_update;
    update.residual = update.mean_action - action_.parameters().S0;
    update.previous_a = action_.parameters().a;
    Updater(update, parameters_.gain, action_.parameters().sigma, status_.iteration);

    action_.set_a(update.updated_a);
    status_.last_update = update;
    ++status_.iteration;
    status_.accumulated_action = 0.0;

    // Swap over the thermalisation steps to adjust for new value of a.
    status_.phase = LLRSolverPhase::Thermalising;
    status_.trajectories_remaining_in_phase = parameters_.rethermalisation_trajectories; 
  }

  /*! @brief Call the thermalise phase of the LLR solver.
   * 
   * Only role here is to switch state to 'Accumulating' when
   * the number of thermalisation steps have passed.
   */
  void thermalise()
  {
    if (status_.trajectories_remaining_in_phase <= 0)
    {
      status_.phase = LLRSolverPhase::Accumulating;
      status_.trajectories_remaining_in_phase = parameters_.trajectories_per_update; 
    }
  }
};

/// \cond DO_NOT_DOCUMENT
NAMESPACE_END(Grid);
/// \endcond

