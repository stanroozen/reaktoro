// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright © 2014-2024 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

// C++ includes
#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <thread>
#include <vector>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/ChemicalProps.hpp>
#include <Reaktoro/Equilibrium/EquilibriumConditions.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSolver.hpp>
#include <Reaktoro/Equilibrium/EquilibriumSweepSolver.hpp>
#include <Reaktoro/Utils/AqueousProps.hpp>
using namespace Reaktoro;

namespace {
auto flattenIndex(Index ix, Index iy, Index ny) -> size_t
{
    return static_cast<size_t>(ix * ny + iy);
}
}

auto EquilibriumSweepResult::size() const -> Index
{
    return static_cast<Index>(results.size());
}

auto EquilibriumSweepResult::succeededCount() const -> Index
{
    Index count = 0;
    for(auto const& result : results)
        count += static_cast<Index>(result.succeeded());
    return count;
}

auto EquilibriumSweepResult::elementMolalityArray(StringOrIndex const& element) const -> ArrayXd
{
    const auto n = static_cast<Index>(states.size());
    ArrayXd out = ArrayXd::Constant(static_cast<Eigen::Index>(n), std::numeric_limits<double>::quiet_NaN());
    for(Index i = 0; i < n; ++i)
    {
        if(!results[static_cast<size_t>(i)].succeeded())
            continue;
        AqueousProps props(states[static_cast<size_t>(i)]);
        out[static_cast<Eigen::Index>(i)] = static_cast<double>(props.elementMolality(element));
    }
    return out;
}

auto EquilibriumSweepGridResult::sizeX() const -> Index
{
    return static_cast<Index>(xvalues.size());
}

auto EquilibriumSweepGridResult::sizeY() const -> Index
{
    return static_cast<Index>(yvalues.size());
}

auto EquilibriumSweepGridResult::size() const -> Index
{
    return static_cast<Index>(results.size());
}

auto EquilibriumSweepGridResult::succeededCount() const -> Index
{
    Index count = 0;
    for(auto const& result : results)
        count += static_cast<Index>(result.succeeded());
    return count;
}

auto EquilibriumSweepGridResult::logActivityGrid(StringOrIndex const& species) const -> ArrayXXd
{
    const auto nx = sizeX();
    const auto ny = sizeY();
    ArrayXXd out = ArrayXXd::Constant(static_cast<Eigen::Index>(nx), static_cast<Eigen::Index>(ny), std::numeric_limits<double>::quiet_NaN());

    if(nx == 0 || ny == 0)
        return out;

    for(Index ix = 0; ix < nx; ++ix)
    {
        for(Index iy = 0; iy < ny; ++iy)
        {
            const auto i = flattenIndex(ix, iy, ny);
            if(!results[i].succeeded())
                continue;

            ChemicalProps props(states[i].system());
            props.update(states[i]);
            out(static_cast<Eigen::Index>(ix), static_cast<Eigen::Index>(iy)) = static_cast<double>(props.speciesActivityLg(species));
        }
    }

    return out;
}

auto EquilibriumSweepGridResult::predominantSpeciesGrid(Strings const& species) const -> ArrayXXd
{
    const auto nx = sizeX();
    const auto ny = sizeY();
    ArrayXXd out = ArrayXXd::Constant(static_cast<Eigen::Index>(nx), static_cast<Eigen::Index>(ny), std::numeric_limits<double>::quiet_NaN());

    if(nx == 0 || ny == 0 || species.empty())
        return out;

    for(Index ix = 0; ix < nx; ++ix)
    {
        for(Index iy = 0; iy < ny; ++iy)
        {
            const auto i = flattenIndex(ix, iy, ny);
            if(!results[i].succeeded())
                continue;

            ChemicalProps props(states[i].system());
            props.update(states[i]);

            auto best_idx = Index(0);
            auto best_lg = -std::numeric_limits<double>::infinity();
            for(Index ispecies = 0; ispecies < static_cast<Index>(species.size()); ++ispecies)
            {
                const auto lg = static_cast<double>(props.speciesActivityLg(species[ispecies]));
                if(std::isnan(lg))
                    continue;
                if(lg > best_lg)
                {
                    best_lg = lg;
                    best_idx = ispecies;
                }
            }

            if(best_lg > -std::numeric_limits<double>::infinity())
                out(static_cast<Eigen::Index>(ix), static_cast<Eigen::Index>(iy)) = static_cast<double>(best_idx);
        }
    }

    return out;
}

auto EquilibriumSweepGridResult::saturationIndexGrid(StringOrIndex const& species) const -> ArrayXXd
{
    const auto nx = sizeX();
    const auto ny = sizeY();
    ArrayXXd out = ArrayXXd::Constant(static_cast<Eigen::Index>(nx), static_cast<Eigen::Index>(ny), std::numeric_limits<double>::quiet_NaN());

    if(nx == 0 || ny == 0)
        return out;

    for(Index ix = 0; ix < nx; ++ix)
    {
        for(Index iy = 0; iy < ny; ++iy)
        {
            const auto i = flattenIndex(ix, iy, ny);
            if(!results[i].succeeded())
                continue;

            AqueousProps props(states[i]);
            out(static_cast<Eigen::Index>(ix), static_cast<Eigen::Index>(iy)) = static_cast<double>(props.saturationIndex(species));
        }
    }

    return out;
}

auto EquilibriumSweepGridResult::elementMolalityGrid(StringOrIndex const& element) const -> ArrayXXd
{
    const auto nx = sizeX();
    const auto ny = sizeY();
    ArrayXXd out = ArrayXXd::Constant(static_cast<Eigen::Index>(nx), static_cast<Eigen::Index>(ny), std::numeric_limits<double>::quiet_NaN());

    if(nx == 0 || ny == 0)
        return out;

    for(Index ix = 0; ix < nx; ++ix)
    {
        for(Index iy = 0; iy < ny; ++iy)
        {
            const auto i = flattenIndex(ix, iy, ny);
            if(!results[i].succeeded())
                continue;

            AqueousProps props(states[i]);
            out(static_cast<Eigen::Index>(ix), static_cast<Eigen::Index>(iy)) = static_cast<double>(props.elementMolality(element));
        }
    }

    return out;
}

EquilibriumSweepSolver::EquilibriumSweepSolver(ChemicalSystem const& system)
    : m_specs(system)
{
}

EquilibriumSweepSolver::EquilibriumSweepSolver(EquilibriumSpecs const& specs)
    : m_specs(specs)
{
}

auto EquilibriumSweepSolver::setOptions(EquilibriumOptions const& options) -> void
{
    m_options = options;
}

auto EquilibriumSweepSolver::setSweepOptions(EquilibriumSweepOptions const& options) -> void
{
    m_sweep_options = options;
}

auto EquilibriumSweepSolver::options() const -> EquilibriumOptions const&
{
    return m_options;
}

auto EquilibriumSweepSolver::sweepOptions() const -> EquilibriumSweepOptions const&
{
    return m_sweep_options;
}

auto EquilibriumSweepSolver::sweepTP(
    ChemicalState const& initial,
    ArrayXdConstRef const& temperatures,
    ArrayXdConstRef const& pressures,
    String const& temperature_unit,
    String const& pressure_unit
) -> EquilibriumSweepResult
{
    const auto n = static_cast<Index>(temperatures.size());
    errorif(n != static_cast<Index>(pressures.size()), "Cannot execute TP sweep because temperature and pressure arrays have different sizes.");

    EquilibriumSweepResult out;
    out.states.assign(static_cast<size_t>(n), initial);
    out.results.resize(static_cast<size_t>(n));

    if(n == 0)
        return out;

    std::vector<double> tvals(static_cast<size_t>(n));
    std::vector<double> pvals(static_cast<size_t>(n));
    for(Index i = 0; i < n; ++i)
    {
        tvals[static_cast<size_t>(i)] = temperatures[i];
        pvals[static_cast<size_t>(i)] = pressures[i];
    }

    const auto num_threads = std::max<Index>(1, m_sweep_options.num_threads);
    if(num_threads <= 1)
    {
        EquilibriumSolver solver(m_specs);
        solver.setOptions(m_options);
        EquilibriumConditions conditions(m_specs);
        ChemicalState state(initial);

        for(Index i = 0; i < n; ++i)
        {
            if(i == 0 || !m_sweep_options.continuation)
                state = initial;

            conditions.temperature(tvals[static_cast<size_t>(i)], temperature_unit);
            conditions.pressure(pvals[static_cast<size_t>(i)], pressure_unit);
            out.results[static_cast<size_t>(i)] = solver.solve(state, conditions);
            out.states[static_cast<size_t>(i)] = state;
        }
        return out;
    }

    const auto workers = static_cast<size_t>(std::min<Index>(num_threads, n));
    std::atomic<Index> next(0);
    std::vector<std::thread> threads;
    threads.reserve(workers);

    auto task = [&]() {
        EquilibriumSolver solver(m_specs);
        solver.setOptions(m_options);
        EquilibriumConditions conditions(m_specs);
        ChemicalState state(initial);

        while(true)
        {
            const auto i = next.fetch_add(1);
            if(i >= n)
                break;

            state = initial;
            conditions.temperature(tvals[static_cast<size_t>(i)], temperature_unit);
            conditions.pressure(pvals[static_cast<size_t>(i)], pressure_unit);
            out.results[static_cast<size_t>(i)] = solver.solve(state, conditions);
            out.states[static_cast<size_t>(i)] = state;
        }
    };

    for(size_t worker = 0; worker < workers; ++worker)
        threads.emplace_back(task);
    for(auto& thread : threads)
        thread.join();

    return out;
}

auto EquilibriumSweepSolver::sweepInput(
    ChemicalState const& initial,
    String const& input,
    ArrayXdConstRef const& values
) -> EquilibriumSweepResult
{
    const auto n = static_cast<Index>(values.size());

    EquilibriumSweepResult out;
    out.states.assign(static_cast<size_t>(n), initial);
    out.results.resize(static_cast<size_t>(n));

    if(n == 0)
        return out;

    std::vector<double> xvals(static_cast<size_t>(n));
    for(Index i = 0; i < n; ++i)
        xvals[static_cast<size_t>(i)] = values[i];

    const auto num_threads = std::max<Index>(1, m_sweep_options.num_threads);
    if(num_threads <= 1)
    {
        EquilibriumSolver solver(m_specs);
        solver.setOptions(m_options);
        EquilibriumConditions conditions(m_specs);
        ChemicalState state(initial);

        for(Index i = 0; i < n; ++i)
        {
            if(i == 0 || !m_sweep_options.continuation)
                state = initial;

            conditions.setInputVariable(input, xvals[static_cast<size_t>(i)]);
            out.results[static_cast<size_t>(i)] = solver.solve(state, conditions);
            out.states[static_cast<size_t>(i)] = state;
        }
        return out;
    }

    const auto workers = static_cast<size_t>(std::min<Index>(num_threads, n));
    std::atomic<Index> next(0);
    std::vector<std::thread> threads;
    threads.reserve(workers);

    auto task = [&]() {
        EquilibriumSolver solver(m_specs);
        solver.setOptions(m_options);
        EquilibriumConditions conditions(m_specs);
        ChemicalState state(initial);

        while(true)
        {
            const auto i = next.fetch_add(1);
            if(i >= n)
                break;

            state = initial;
            conditions.setInputVariable(input, xvals[static_cast<size_t>(i)]);
            out.results[static_cast<size_t>(i)] = solver.solve(state, conditions);
            out.states[static_cast<size_t>(i)] = state;
        }
    };

    for(size_t worker = 0; worker < workers; ++worker)
        threads.emplace_back(task);
    for(auto& thread : threads)
        thread.join();

    return out;
}

auto EquilibriumSweepSolver::sweepPH(
    ChemicalState const& initial,
    ArrayXdConstRef const& values
) -> EquilibriumSweepResult
{
    const auto n = static_cast<Index>(values.size());

    EquilibriumSweepResult out;
    out.states.assign(static_cast<size_t>(n), initial);
    out.results.resize(static_cast<size_t>(n));

    if(n == 0)
        return out;

    std::vector<double> xvals(static_cast<size_t>(n));
    for(Index i = 0; i < n; ++i)
        xvals[static_cast<size_t>(i)] = values[i];

    const auto num_threads = std::max<Index>(1, m_sweep_options.num_threads);
    if(num_threads <= 1)
    {
        EquilibriumSolver solver(m_specs);
        solver.setOptions(m_options);
        EquilibriumConditions conditions(m_specs);
        ChemicalState state(initial);

        for(Index i = 0; i < n; ++i)
        {
            if(i == 0 || !m_sweep_options.continuation)
                state = initial;

            conditions.pH(xvals[static_cast<size_t>(i)]);
            out.results[static_cast<size_t>(i)] = solver.solve(state, conditions);
            out.states[static_cast<size_t>(i)] = state;
        }
        return out;
    }

    const auto workers = static_cast<size_t>(std::min<Index>(num_threads, n));
    std::atomic<Index> next(0);
    std::vector<std::thread> threads;
    threads.reserve(workers);

    auto task = [&]() {
        EquilibriumSolver solver(m_specs);
        solver.setOptions(m_options);
        EquilibriumConditions conditions(m_specs);
        ChemicalState state(initial);

        while(true)
        {
            const auto i = next.fetch_add(1);
            if(i >= n)
                break;

            state = initial;
            conditions.pH(xvals[static_cast<size_t>(i)]);
            out.results[static_cast<size_t>(i)] = solver.solve(state, conditions);
            out.states[static_cast<size_t>(i)] = state;
        }
    };

    for(size_t worker = 0; worker < workers; ++worker)
        threads.emplace_back(task);
    for(auto& thread : threads)
        thread.join();

    return out;
}

auto EquilibriumSweepSolver::sweepEh(
    ChemicalState const& initial,
    ArrayXdConstRef const& values,
    String const& unit
) -> EquilibriumSweepResult
{
    const auto n = static_cast<Index>(values.size());

    EquilibriumSweepResult out;
    out.states.assign(static_cast<size_t>(n), initial);
    out.results.resize(static_cast<size_t>(n));

    if(n == 0)
        return out;

    std::vector<double> xvals(static_cast<size_t>(n));
    for(Index i = 0; i < n; ++i)
        xvals[static_cast<size_t>(i)] = values[i];

    const auto num_threads = std::max<Index>(1, m_sweep_options.num_threads);
    if(num_threads <= 1)
    {
        EquilibriumSolver solver(m_specs);
        solver.setOptions(m_options);
        EquilibriumConditions conditions(m_specs);
        ChemicalState state(initial);

        for(Index i = 0; i < n; ++i)
        {
            if(i == 0 || !m_sweep_options.continuation)
                state = initial;

            conditions.Eh(xvals[static_cast<size_t>(i)], unit);
            out.results[static_cast<size_t>(i)] = solver.solve(state, conditions);
            out.states[static_cast<size_t>(i)] = state;
        }
        return out;
    }

    const auto workers = static_cast<size_t>(std::min<Index>(num_threads, n));
    std::atomic<Index> next(0);
    std::vector<std::thread> threads;
    threads.reserve(workers);

    auto task = [&]() {
        EquilibriumSolver solver(m_specs);
        solver.setOptions(m_options);
        EquilibriumConditions conditions(m_specs);
        ChemicalState state(initial);

        while(true)
        {
            const auto i = next.fetch_add(1);
            if(i >= n)
                break;

            state = initial;
            conditions.Eh(xvals[static_cast<size_t>(i)], unit);
            out.results[static_cast<size_t>(i)] = solver.solve(state, conditions);
            out.states[static_cast<size_t>(i)] = state;
        }
    };

    for(size_t worker = 0; worker < workers; ++worker)
        threads.emplace_back(task);
    for(auto& thread : threads)
        thread.join();

    return out;
}

auto EquilibriumSweepSolver::sweepPHEhGrid(
    ChemicalState const& initial,
    ArrayXdConstRef const& pH_values,
    ArrayXdConstRef const& Eh_values,
    String const& Eh_unit
) -> EquilibriumSweepGridResult
{
    const auto nx = static_cast<Index>(pH_values.size());
    const auto ny = static_cast<Index>(Eh_values.size());
    const auto n = nx * ny;

    EquilibriumSweepGridResult out;
    out.xvalues = pH_values;
    out.yvalues = Eh_values;
    out.states.assign(static_cast<size_t>(n), initial);
    out.results.resize(static_cast<size_t>(n));

    if(n == 0)
        return out;

    std::vector<double> xvals(static_cast<size_t>(nx));
    std::vector<double> yvals(static_cast<size_t>(ny));
    for(Index ix = 0; ix < nx; ++ix)
        xvals[static_cast<size_t>(ix)] = pH_values[ix];
    for(Index iy = 0; iy < ny; ++iy)
        yvals[static_cast<size_t>(iy)] = Eh_values[iy];

    const auto num_threads = std::max<Index>(1, m_sweep_options.num_threads);
    if(num_threads <= 1)
    {
        EquilibriumSolver solver(m_specs);
        solver.setOptions(m_options);
        EquilibriumConditions conditions(m_specs);
        ChemicalState state(initial);

        for(Index ix = 0; ix < nx; ++ix)
        {
            for(Index iy = 0; iy < ny; ++iy)
            {
                if((ix == 0 && iy == 0) || !m_sweep_options.continuation)
                    state = initial;

                conditions.pH(xvals[static_cast<size_t>(ix)]);
                conditions.Eh(yvals[static_cast<size_t>(iy)], Eh_unit);
                const auto i = flattenIndex(ix, iy, ny);
                out.results[i] = solver.solve(state, conditions);
                out.states[i] = state;
            }
        }
        return out;
    }

    const auto workers = static_cast<size_t>(std::min<Index>(num_threads, n));
    std::atomic<Index> next(0);
    std::vector<std::thread> threads;
    threads.reserve(workers);

    auto task = [&]() {
        EquilibriumSolver solver(m_specs);
        solver.setOptions(m_options);
        EquilibriumConditions conditions(m_specs);
        ChemicalState state(initial);

        while(true)
        {
            const auto i = next.fetch_add(1);
            if(i >= n)
                break;

            const auto ix = i / ny;
            const auto iy = i % ny;

            state = initial;
            conditions.pH(xvals[static_cast<size_t>(ix)]);
            conditions.Eh(yvals[static_cast<size_t>(iy)], Eh_unit);
            out.results[static_cast<size_t>(i)] = solver.solve(state, conditions);
            out.states[static_cast<size_t>(i)] = state;
        }
    };

    for(size_t worker = 0; worker < workers; ++worker)
        threads.emplace_back(task);
    for(auto& thread : threads)
        thread.join();

    return out;
}

auto EquilibriumSweepSolver::sweepLgActivityGrid(
    ChemicalState const& initial,
    String const& speciesX,
    ArrayXdConstRef const& lgaX_values,
    String const& speciesY,
    ArrayXdConstRef const& lgaY_values
) -> EquilibriumSweepGridResult
{
    const auto nx = static_cast<Index>(lgaX_values.size());
    const auto ny = static_cast<Index>(lgaY_values.size());
    const auto n = nx * ny;

    EquilibriumSweepGridResult out;
    out.xvalues = lgaX_values;
    out.yvalues = lgaY_values;
    out.states.assign(static_cast<size_t>(n), initial);
    out.results.resize(static_cast<size_t>(n));

    if(n == 0)
        return out;

    std::vector<double> xvals(static_cast<size_t>(nx));
    std::vector<double> yvals(static_cast<size_t>(ny));
    for(Index ix = 0; ix < nx; ++ix)
        xvals[static_cast<size_t>(ix)] = lgaX_values[ix];
    for(Index iy = 0; iy < ny; ++iy)
        yvals[static_cast<size_t>(iy)] = lgaY_values[iy];

    const auto num_threads = std::max<Index>(1, m_sweep_options.num_threads);
    if(num_threads <= 1)
    {
        EquilibriumSolver solver(m_specs);
        solver.setOptions(m_options);
        EquilibriumConditions conditions(m_specs);
        ChemicalState state(initial);

        for(Index ix = 0; ix < nx; ++ix)
        {
            for(Index iy = 0; iy < ny; ++iy)
            {
                if((ix == 0 && iy == 0) || !m_sweep_options.continuation)
                    state = initial;

                conditions.lgActivity(speciesX, xvals[static_cast<size_t>(ix)]);
                conditions.lgActivity(speciesY, yvals[static_cast<size_t>(iy)]);
                const auto i = flattenIndex(ix, iy, ny);
                out.results[i] = solver.solve(state, conditions);
                out.states[i] = state;
            }
        }
        return out;
    }

    const auto workers = static_cast<size_t>(std::min<Index>(num_threads, n));
    std::atomic<Index> next(0);
    std::vector<std::thread> threads;
    threads.reserve(workers);

    auto task = [&]() {
        EquilibriumSolver solver(m_specs);
        solver.setOptions(m_options);
        EquilibriumConditions conditions(m_specs);
        ChemicalState state(initial);

        while(true)
        {
            const auto i = next.fetch_add(1);
            if(i >= n)
                break;

            const auto ix = i / ny;
            const auto iy = i % ny;

            state = initial;
            conditions.lgActivity(speciesX, xvals[static_cast<size_t>(ix)]);
            conditions.lgActivity(speciesY, yvals[static_cast<size_t>(iy)]);
            out.results[static_cast<size_t>(i)] = solver.solve(state, conditions);
            out.states[static_cast<size_t>(i)] = state;
        }
    };

    for(size_t worker = 0; worker < workers; ++worker)
        threads.emplace_back(task);
    for(auto& thread : threads)
        thread.join();

    return out;
}

auto EquilibriumSweepSolver::sweepTPGrid(
    ChemicalState const& initial,
    ArrayXdConstRef const& temperature_values,
    String const& temperature_unit,
    ArrayXdConstRef const& pressure_values,
    String const& pressure_unit
) -> EquilibriumSweepGridResult
{
    const auto nx = static_cast<Index>(temperature_values.size());
    const auto ny = static_cast<Index>(pressure_values.size());
    const auto n = nx * ny;

    EquilibriumSweepGridResult out;
    out.xvalues = temperature_values;
    out.yvalues = pressure_values;
    out.states.assign(static_cast<size_t>(n), initial);
    out.results.resize(static_cast<size_t>(n));

    if(n == 0)
        return out;

    std::vector<double> Tvals(static_cast<size_t>(nx));
    std::vector<double> Pvals(static_cast<size_t>(ny));
    for(Index ix = 0; ix < nx; ++ix)
        Tvals[static_cast<size_t>(ix)] = temperature_values[ix];
    for(Index iy = 0; iy < ny; ++iy)
        Pvals[static_cast<size_t>(iy)] = pressure_values[iy];

    const auto num_threads = std::max<Index>(1, m_sweep_options.num_threads);
    if(num_threads <= 1)
    {
        EquilibriumSolver solver(m_specs);
        solver.setOptions(m_options);
        EquilibriumConditions conditions(m_specs);
        ChemicalState state(initial);

        for(Index ix = 0; ix < nx; ++ix)
        {
            for(Index iy = 0; iy < ny; ++iy)
            {
                if((ix == 0 && iy == 0) || !m_sweep_options.continuation)
                    state = initial;

                conditions.temperature(Tvals[static_cast<size_t>(ix)], temperature_unit);
                conditions.pressure(Pvals[static_cast<size_t>(iy)], pressure_unit);
                const auto i = flattenIndex(ix, iy, ny);
                out.results[i] = solver.solve(state, conditions);
                out.states[i] = state;
            }
        }
        return out;
    }

    const auto workers = static_cast<size_t>(std::min<Index>(num_threads, n));
    std::atomic<Index> next(0);
    std::vector<std::thread> threads;
    threads.reserve(workers);

    auto task = [&]() {
        EquilibriumSolver solver(m_specs);
        solver.setOptions(m_options);
        EquilibriumConditions conditions(m_specs);
        ChemicalState state(initial);

        while(true)
        {
            const auto i = next.fetch_add(1);
            if(i >= n)
                break;

            const auto ix = i / ny;
            const auto iy = i % ny;

            state = initial;
            conditions.temperature(Tvals[static_cast<size_t>(ix)], temperature_unit);
            conditions.pressure(Pvals[static_cast<size_t>(iy)], pressure_unit);
            out.results[static_cast<size_t>(i)] = solver.solve(state, conditions);
            out.states[static_cast<size_t>(i)] = state;
        }
    };

    for(size_t worker = 0; worker < workers; ++worker)
        threads.emplace_back(task);
    for(auto& thread : threads)
        thread.join();

    return out;
}

auto EquilibriumSweepSolver::sweepLogfO2pHGrid(
    ChemicalState const& initial,
    ArrayXdConstRef const& logfO2_values,
    String const& fug_unit,
    ArrayXdConstRef const& pH_values
) -> EquilibriumSweepGridResult
{
    const auto nx = static_cast<Index>(logfO2_values.size());
    const auto ny = static_cast<Index>(pH_values.size());
    const auto n = nx * ny;

    EquilibriumSweepGridResult out;
    out.xvalues = logfO2_values;
    out.yvalues = pH_values;
    out.states.assign(static_cast<size_t>(n), initial);
    out.results.resize(static_cast<size_t>(n));

    if(n == 0)
        return out;

    std::vector<double> lfO2(static_cast<size_t>(nx));
    std::vector<double> pHvals(static_cast<size_t>(ny));
    for(Index ix = 0; ix < nx; ++ix)
        lfO2[static_cast<size_t>(ix)] = logfO2_values[ix];
    for(Index iy = 0; iy < ny; ++iy)
        pHvals[static_cast<size_t>(iy)] = pH_values[iy];

    const auto num_threads = std::max<Index>(1, m_sweep_options.num_threads);
    if(num_threads <= 1)
    {
        EquilibriumSolver solver(m_specs);
        solver.setOptions(m_options);
        EquilibriumConditions conditions(m_specs);
        ChemicalState state(initial);

        for(Index ix = 0; ix < nx; ++ix)
        {
            for(Index iy = 0; iy < ny; ++iy)
            {
                if((ix == 0 && iy == 0) || !m_sweep_options.continuation)
                    state = initial;

                const double fO2_val = std::pow(10.0, lfO2[static_cast<size_t>(ix)]);
                conditions.fugacity("O2", fO2_val, fug_unit);
                conditions.pH(pHvals[static_cast<size_t>(iy)]);
                const auto i = flattenIndex(ix, iy, ny);
                out.results[i] = solver.solve(state, conditions);
                out.states[i] = state;
            }
        }
        return out;
    }

    const auto workers = static_cast<size_t>(std::min<Index>(num_threads, n));
    std::atomic<Index> next(0);
    std::vector<std::thread> threads;
    threads.reserve(workers);

    auto task = [&]() {
        EquilibriumSolver solver(m_specs);
        solver.setOptions(m_options);
        EquilibriumConditions conditions(m_specs);
        ChemicalState state(initial);

        while(true)
        {
            const auto i = next.fetch_add(1);
            if(i >= n)
                break;

            const auto ix = i / ny;
            const auto iy = i % ny;

            state = initial;
            const double fO2_val = std::pow(10.0, lfO2[static_cast<size_t>(ix)]);
            conditions.fugacity("O2", fO2_val, fug_unit);
            conditions.pH(pHvals[static_cast<size_t>(iy)]);
            out.results[static_cast<size_t>(i)] = solver.solve(state, conditions);
            out.states[static_cast<size_t>(i)] = state;
        }
    };

    for(size_t worker = 0; worker < workers; ++worker)
        threads.emplace_back(task);
    for(auto& thread : threads)
        thread.join();

    return out;
}
