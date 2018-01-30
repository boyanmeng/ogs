/**
 * \author Norihiro Watanabe
 * \date   2012-08-03
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <cmath>
#include <vector>

#include "NumLib/TimeStepping/TimeStep.h"

namespace NumLib
{
/**
 * \brief Interface of time stepping algorithms
 *
 */
class TimeStepAlgorithm
{
public:
    TimeStepAlgorithm(const double t0, const double t_end)
        : _t_initial(t0), _t_end(t_end), _ts_prev(t0), _ts_current(t0)
    {
    }

    TimeStepAlgorithm(const double t0, const double t_end, const double dt)
        : _t_initial(t0),
          _t_end(t_end),
          _ts_prev(t0),
          _ts_current(t0),
          _dt_vector(static_cast<std::size_t>(std::ceil((t_end - t0) / dt)), dt)
    {
    }

    TimeStepAlgorithm(const double t0, const double t_end,
                      const std::vector<double>& all_step_sizes)
        : _t_initial(t0),
          _t_end(t_end),
          _ts_prev(t0),
          _ts_current(t0),
          _dt_vector(all_step_sizes)
    {
    }

    virtual ~TimeStepAlgorithm() = default;

    /// return the beginning of time steps
    double begin() const { return _t_initial; }
    /// return the end of time steps
    double end() const { return _t_end; }
    /// return current time step
    const TimeStep getTimeStep() const { return _ts_current; }
    /// reset the current step size from the previous time
    void resetCurrentTimeStep(const double dt)
    {
        _ts_current = _ts_prev;
        _ts_current += dt;
    }

    /// move to the next time step
    /// \param solution_error Solution error between two successive time steps.
    /// \return true if the next step exists
    virtual bool next(const double solution_error) = 0;

    /// return if current time step is accepted or not
    virtual bool accepted() const = 0;

    /// Set the status of the step. A boolean parameter is needed to indicated
    /// whether the step is accepted or not.
    virtual void setAcceptedOrNot(const bool /*accepted*/){};

    /// return a history of time step sizes
    const std::vector<double>& getTimeStepSizeHistory() const
    {
        return _dt_vector;
    }

    /// Get a flag to indicate whether this algorithm needs to compute
    /// solution error. The default return value is false.
    virtual bool isSolutionErrorComputationNeeded() { return false; }
protected:
    /// initial time
    const double _t_initial;
    /// end time
    const double _t_end;

    /// previous time step information
    TimeStep _ts_prev;
    /// current time step information
    TimeStep _ts_current;

    /// a vector of time step sizes
    std::vector<double> _dt_vector;
};

}  // NumLib
