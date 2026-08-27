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
#include <utility>
#include <vector>

// nlohmann includes
#include <nlohmann/json.hpp>

// Reaktoro includes
#include <Reaktoro/Common/Exception.hpp>
#include <Reaktoro/Core/AggregateState.hpp>
#include <Reaktoro/Core/Database.hpp>
#include <Reaktoro/Equilibrium/EquilibriumBenchmarkUtils.hpp>
using namespace Reaktoro;

namespace {

auto nanValue() -> double
{
    return std::numeric_limits<double>::quiet_NaN();
}

auto clamp(double value, double lo, double hi) -> double
{
    return std::max(lo, std::min(hi, value));
}

auto quantileSorted(std::vector<double> const& sorted, double q) -> double
{
    if(sorted.empty())
        return nanValue();

    const auto qq = clamp(q, 0.0, 1.0);
    const auto n = sorted.size();
    if(n == 1)
        return sorted[0];

    const double h = qq * static_cast<double>(n - 1);
    const auto i0 = static_cast<size_t>(std::floor(h));
    const auto i1 = static_cast<size_t>(std::ceil(h));
    if(i0 == i1)
        return sorted[i0];

    const double w = h - static_cast<double>(i0);
    return sorted[i0] * (1.0 - w) + sorted[i1] * w;
}

auto collectFinite(ArrayXdConstRef const& x, ArrayXdConstRef const& y) -> std::vector<std::pair<double, double>>
{
    std::vector<std::pair<double, double>> xy;
    const auto n = static_cast<Index>(x.size());
    xy.reserve(static_cast<size_t>(n));
    for(Index i = 0; i < n; ++i)
    {
        const auto xv = static_cast<double>(x[i]);
        const auto yv = static_cast<double>(y[i]);
        if(std::isfinite(xv) && std::isfinite(yv))
            xy.emplace_back(xv, yv);
    }
    return xy;
}

} // namespace

auto Reaktoro::perturbMineralDatabaseJSON(
    String const& base_json,
    Strings const& entities,
    ArrayXdConstRef const& shifts_j_per_mol
) -> String
{
    const auto n = static_cast<Index>(entities.size());
    errorif(shifts_j_per_mol.size() != n, "Cannot perturb mineral JSON because entity and shift sizes differ.");

    auto sampled = nlohmann::json::parse(base_json);
    auto& species = sampled["Species"];

    for(Index i = 0; i < n; ++i)
    {
        const auto& code = entities[static_cast<size_t>(i)];
        if(!species.contains(code))
            continue;
        auto& sp = species[code];
        if(!sp.contains("StandardThermoModel"))
            continue;
        auto& stdm = sp["StandardThermoModel"];
        if(!stdm.contains("HollandPowell"))
            continue;
        auto& hp = stdm["HollandPowell"];
        if(!hp.contains("Gf"))
            continue;

        const auto g0 = hp["Gf"].get<double>();
        hp["Gf"] = g0 + static_cast<double>(shifts_j_per_mol[i]);
    }

    return sampled.dump();
}

auto Reaktoro::perturbMineralDatabase(
    String const& base_json,
    Strings const& entities,
    ArrayXdConstRef const& shifts_j_per_mol
) -> Database
{
    const auto sampled_json = perturbMineralDatabaseJSON(base_json, entities, shifts_j_per_mol);
    return Database::fromStringJSON(sampled_json);
}

auto Reaktoro::perturbMineralDatabases(
    String const& base_json,
    Strings const& entities,
    MatrixXdConstRef const& shifts_j_per_mol_samples,
    Index num_threads
) -> Vec<Database>
{
    const auto nsamples = static_cast<Index>(shifts_j_per_mol_samples.rows());
    const auto nentities = static_cast<Index>(entities.size());
    errorif(static_cast<Index>(shifts_j_per_mol_samples.cols()) != nentities,
        "Cannot perturb mineral databases because matrix columns do not match number of entities.");

    Vec<Database> databases(static_cast<size_t>(nsamples));
    if(nsamples == 0)
        return databases;

    const auto workers = std::max<Index>(1, std::min(num_threads, nsamples));
    if(workers == 1)
    {
        for(Index i = 0; i < nsamples; ++i)
        {
            const auto shifts = shifts_j_per_mol_samples.row(i).transpose().array();
            databases[static_cast<size_t>(i)] = perturbMineralDatabase(base_json, entities, shifts);
        }
        return databases;
    }

    std::atomic<Index> next(0);
    std::vector<std::thread> threads;
    threads.reserve(static_cast<size_t>(workers));

    auto task = [&]() {
        while(true)
        {
            const auto i = next.fetch_add(1);
            if(i >= nsamples)
                break;
            const auto shifts = shifts_j_per_mol_samples.row(i).transpose().array();
            databases[static_cast<size_t>(i)] = perturbMineralDatabase(base_json, entities, shifts);
        }
    };

    for(Index w = 0; w < workers; ++w)
        threads.emplace_back(task);
    for(auto& t : threads)
        t.join();

    return databases;
}

auto Reaktoro::aqueousSpeciesNamesWithAllowedElements(
    Database const& database,
    Strings const& allowed_elements,
    Strings const& excluded_species
) -> Strings
{
    Set<String> allowed;
    for(auto const& symbol : allowed_elements)
        allowed.insert(symbol);

    Set<String> excluded;
    for(auto const& name : excluded_species)
        excluded.insert(name);

    Strings names;
    const auto aqueous = database.speciesWithAggregateState(AggregateState::Aqueous);
    names.reserve(static_cast<size_t>(aqueous.size()));

    for(auto const& species : aqueous)
    {
        const auto name = species.name();
        if(excluded.find(name) != excluded.end())
            continue;

        bool ok = true;
        for(auto const& elem : species.formula().elements())
        {
            if(allowed.find(elem.first) == allowed.end())
            {
                ok = false;
                break;
            }
        }

        if(ok)
            names.push_back(name);
    }

    return names;
}

auto Reaktoro::elementStoichiometryTerms(
    Database const& database,
    Strings const& species_names,
    Strings const& target_elements
) -> Vec<ElementStoichiometryTerm>
{
    Vec<ElementStoichiometryTerm> out;
    for(auto const& species_name : species_names)
    {
        const auto species = database.species(species_name);
        const auto formula = species.formula();

        for(auto const& element : target_elements)
        {
            const auto coeff = formula.coefficient(element);
            if(coeff == 0.0)
                continue;
            out.push_back({element, species_name, coeff});
        }
    }
    return out;
}

auto Reaktoro::interpolateCurveValue(
    ArrayXdConstRef const& x,
    ArrayXdConstRef const& y,
    double x_query,
    double atol
) -> double
{
    errorif(x.size() != y.size(), "Cannot interpolate curve value because x and y arrays have different sizes.");

    auto xy = collectFinite(x, y);
    if(xy.empty())
        return nanValue();

    for(auto const& p : xy)
    {
        if(std::abs(p.first - x_query) <= atol)
            return p.second;
    }

    if(xy.size() == 1)
        return nanValue();

    std::sort(xy.begin(), xy.end(), [](auto const& a, auto const& b) { return a.first < b.first; });

    if(x_query < xy.front().first || x_query > xy.back().first)
        return nanValue();

    auto it = std::lower_bound(
        xy.begin(), xy.end(), x_query,
        [](auto const& a, double value) { return a.first < value; }
    );

    if(it == xy.begin())
        return it->second;
    if(it == xy.end())
        return xy.back().second;
    if(std::abs(it->first - x_query) <= atol)
        return it->second;

    auto const& right = *it;
    auto const& left = *(it - 1);
    const auto dx = right.first - left.first;
    if(std::abs(dx) <= std::numeric_limits<double>::epsilon())
        return left.second;

    const auto w = (x_query - left.first) / dx;
    return left.second + w * (right.second - left.second);
}

auto Reaktoro::computeResidualsInterpolated(
    ArrayXdConstRef const& curve_x,
    ArrayXdConstRef const& curve_y,
    ArrayXdConstRef const& query_x,
    ArrayXdConstRef const& query_y,
    double atol
) -> InterpolationResidualResult
{
    errorif(query_x.size() != query_y.size(), "Cannot compute residuals because query arrays have different sizes.");

    const auto n = static_cast<Index>(query_x.size());
    InterpolationResidualResult out;
    out.calculated = ArrayXd::Constant(static_cast<Eigen::Index>(n), nanValue());
    out.residuals = ArrayXd::Constant(static_cast<Eigen::Index>(n), nanValue());

    for(Index i = 0; i < n; ++i)
    {
        const auto xq = static_cast<double>(query_x[i]);
        const auto yq = static_cast<double>(query_y[i]);
        if(!std::isfinite(xq) || !std::isfinite(yq))
            continue;

        const auto yc = interpolateCurveValue(curve_x, curve_y, xq, atol);
        if(!std::isfinite(yc))
            continue;

        out.calculated[i] = yc;
        out.residuals[i] = yc - yq;
    }

    return out;
}

auto Reaktoro::computeUncertaintyBand(
    ArrayXXdConstRef const& samples,
    double ci_percent
) -> UncertaintyBandResult
{
    const auto nsamples = static_cast<Index>(samples.rows());
    const auto npoints = static_cast<Index>(samples.cols());

    UncertaintyBandResult out;
    out.lower = ArrayXd::Constant(static_cast<Eigen::Index>(npoints), nanValue());
    out.median = ArrayXd::Constant(static_cast<Eigen::Index>(npoints), nanValue());
    out.upper = ArrayXd::Constant(static_cast<Eigen::Index>(npoints), nanValue());

    if(nsamples == 0 || npoints == 0)
        return out;

    const auto ci = clamp(ci_percent, 0.0, 100.0);
    const auto alpha = 0.5 * (1.0 - ci / 100.0);
    const auto qlo = alpha;
    const auto qmed = 0.5;
    const auto qhi = 1.0 - alpha;

    std::vector<double> column;
    column.reserve(static_cast<size_t>(nsamples));

    for(Index j = 0; j < npoints; ++j)
    {
        column.clear();
        for(Index i = 0; i < nsamples; ++i)
        {
            const auto v = static_cast<double>(samples(i, j));
            if(std::isfinite(v))
                column.push_back(v);
        }

        if(column.empty())
            continue;

        std::sort(column.begin(), column.end());
        out.lower[j] = quantileSorted(column, qlo);
        out.median[j] = quantileSorted(column, qmed);
        out.upper[j] = quantileSorted(column, qhi);
    }

    return out;
}
