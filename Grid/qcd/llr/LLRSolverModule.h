/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: LLRSolverModule.h

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
#include <string>

#include <Grid/qcd/llr/LLRSolver.h>
#include <Grid/qcd/observables/hmc_observable.h>
#include <Grid/qcd/modules/ObservableModules.h>

NAMESPACE_BEGIN(Grid);

inline void LogLLRSolverUpdate(const LLRSolverUpdate &update,
			       const std::string &context)
{
  std::cout << GridLogMessage << "[LLRSolver" << context << "] trajectory: " << update.trajectory
            << ", iteration: " << update.iteration
            << ", mean action: " << update.mean_action
            << ", residual: " << update.residual
            << ", a: " << update.previous_a << " -> " << update.updated_a << std::endl;
}

template <class LLRSolver>
class LLRActionLogger : public HmcObservable<typename LLRSolver::Field> {

public:
  typedef typename LLRSolver::Field Field;
  
  LLRActionLogger(LLRSolver &solver) : solver_(solver) {}

  virtual void TrajectoryComplete(int trajectory, Field &U,
                                  GridSerialRNG &, GridParallelRNG &)
  {
    int previous_iteration = solver_.status().iteration;
    solver_.record_configuration(trajectory, U);
    report_update(previous_iteration);
  }

  virtual void TrajectoryComplete(int trajectory,
                                  ConfigurationBase<Field> &configuration,
                                  GridSerialRNG &, GridParallelRNG &,
                                  bool)
  {
    int previous_iteration = solver_.status().iteration;
    solver_.record_configuration(trajectory, configuration);
    report_update(previous_iteration);
  }

  std::vector<std::pair<int, RealD> > getMeanActionVector()
  {
    std::vector<std::pair<int, RealD> > out;
    out = mean_action_;
    return out;
  }

  std::vector<std::pair<int, RealD> > getAVector()
  {
    std::vector<std::pair<int, RealD> > out;
    out = a_;
    return out;
  }


private:
  LLRSolver &solver_;
  std::vector<std::pair<int, RealD> > mean_action_;
  std::vector<std::pair<int, RealD> > a_;

  void report_update(int previous_iteration)
  {
    LLRSolverStatus status = solver_.status();
    if (status.iteration == previous_iteration)
    {
      return;
    }

    LogLLRSolverUpdate(status.last_update, "");
    mean_action_.push_back(std::pair<int, RealD>(status.last_update.trajectory, status.last_update.mean_action));
    a_.push_back(std::pair<int, RealD>(status.last_update.trajectory, status.last_update.previous_a));
    std::cout << GridLogMessage << "[Action Parameters]" << std::endl;
    std::cout << solver_.action_.LogParameters();
  }
};



template <class LLRSolver>
class LLRSolverModule
  : public ObservableModule<LLRActionLogger<LLRSolver>, NoParameters> {
public:
  typedef typename LLRSolver::Field Field;

  typedef ObservableModule<LLRActionLogger<LLRSolver>, NoParameters> ObsBase;
  using ObsBase::ObsBase;

  LLRSolverModule(LLRSolver &solver): ObsBase(NoParameters()) {
    this->ObservablePtr.reset(new LLRActionLogger<LLRSolver>(solver));
  }
  
  virtual void initialize(){}
};

NAMESPACE_END(Grid);

