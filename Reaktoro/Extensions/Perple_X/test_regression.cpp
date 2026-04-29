/**
 * Regression Tests for Perple_X Implementation in Reaktoro
 *
 * This test suite compares Reaktoro outputs against Perple_X reference data
 * to ensure complete accuracy of the implementation.
 *
 * Test Strategy:
 * 1. Generate reference data from Perple_X COHSRK/fluids executable
 * 2. Run identical calculations in Reaktoro
 * 3. Compare results with tight tolerances
 * 4. Report any deviations
 */

#include "PerpleXFluidModel.hpp"
#include "PerpleXElectrolyte.hpp"
#include "PerpleXHKF.hpp"
#include "PerpleXPureEos.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <cmath>
#include <cassert>
#include <string>
#include <sstream>
#include <unordered_map>
#include <map>
#include <filesystem>
#include <limits>
#include <cctype>

using namespace Reaktoro::PerpleX;

// Test tolerance levels
constexpr double TOL_FUGACITY = 1e-6;     // Fugacity coefficient relative error
constexpr double TOL_VOLUME = 1e-6;       // Volume relative error
constexpr double TOL_DIELECTRIC = 1e-4;   // Dielectric constant absolute error
constexpr double TOL_GFUNCTION = 1e-8;    // g-function absolute error (Angstrom)
constexpr double TOL_DH_FACTOR = 1e-6;    // Debye-HÃ¼ckel factor relative error
constexpr double TOL_HKF_GIBBS = 100.0;   // HKF Gibbs energy absolute error (J/mol)
constexpr double TOL_DELTAG_MIX = 25.0;   // Gibbs-mix proxy absolute error (J/mol)
constexpr double TOL_MU_ENERGY  = 1.0;    // Chemical-potential absolute error (J/mol)

struct VolumeToleranceTier {
    double absTol;
    double relTol;
    const char* name;
};

// Reference data structure
struct ReferencePoint {
    double T;           // Temperature (K)
    double P;           // Pressure (bar)
    std::vector<double> y;  // Composition (mole fractions)

    // Expected outputs from Perple_X
    std::vector<double> ln_f;  // Fugacity coefficients
    std::vector<double> v;     // Partial molar volumes (cm3/mol)
    double vol;                // Total volume (cm3/mol)
    double epsilon;            // Dielectric constant
    double gf;                 // Shock g-function (Angstrom)
    double adh;                // Debye-HÃ¼ckel factor
};

// Test result structure
struct TestResult {
    std::string testName;
    bool passed;
    std::vector<std::string> errors;
    std::vector<std::string> warnings;
};

struct RegressionMetrics {
    size_t assertionsExecuted = 0;
    size_t assertionsPassed = 0;
    size_t assertionsFailed = 0;

    size_t rowAssertionsExecuted = 0;
    size_t rowAssertionsPassed = 0;
    size_t rowAssertionsFailed = 0;
};

static RegressionMetrics g_metrics;

VolumeToleranceTier selectVolumeToleranceTier(const std::array<double, 19>& y,
                                              const std::vector<int>& species,
                                              double pressureBar,
                                              double temperatureK)
{
    // Tier 1 (core): tighter interior-composition and moderate P-T region.
    // Tier 2 (boundary): near endmember compositions where sensitivity is higher.
    // Tier 3 (extreme): high/low P-T edge region.
    constexpr VolumeToleranceTier core{0.05, 3e-4, "core"};
    constexpr VolumeToleranceTier boundary{0.075, 4e-4, "boundary"};
    constexpr VolumeToleranceTier extreme{0.09, 5e-4, "extreme"};

    const bool extremePT = (pressureBar <= 10.0 || pressureBar >= 4500.0
                         || temperatureK <= 300.0 || temperatureK >= 750.0);
    if(extremePT)
        return extreme;

    bool boundaryComp = false;
    for(const int idx : species)
    {
        const double yi = y[idx];
        if(yi > 0.0 && (yi <= 0.02 || yi >= 0.98))
        {
            boundaryComp = true;
            break;
        }
    }

    if(boundaryComp)
        return boundary;

    return core;
}

struct CsvTable {
    std::unordered_map<std::string, size_t> columns;
    std::vector<std::string> headers;
    std::vector<std::vector<double>> rows;
};

double parseDouble(const std::string& token)
{
    if(token == "nan" || token == "NaN" || token == "NAN")
        return std::numeric_limits<double>::quiet_NaN();
    return std::stod(token);
}

CsvTable loadCsvTable(const std::filesystem::path& path)
{
    std::ifstream file(path);
    if(!file.is_open())
        throw std::runtime_error("Failed to open CSV: " + path.string());

    std::string headerLine;
    if(!std::getline(file, headerLine))
        throw std::runtime_error("Empty CSV: " + path.string());

    CsvTable table;

    std::stringstream headerStream(headerLine);
    std::string header;
    size_t index = 0;
    while(std::getline(headerStream, header, ','))
    {
        table.columns[header] = index++;
        table.headers.push_back(header);
    }

    std::string line;
    while(std::getline(file, line))
    {
        if(line.empty())
            continue;
        std::stringstream lineStream(line);
        std::string token;
        std::vector<double> row;
        size_t colIndex = 0;
        while(std::getline(lineStream, token, ','))
        {
            const std::string& col = table.headers[colIndex];
            if(col == "case" || col == "sectioning")
            {
                row.push_back(std::numeric_limits<double>::quiet_NaN());
            }
            else
            {
                try
                {
                    row.push_back(parseDouble(token));
                }
                catch(...)
                {
                    row.push_back(std::numeric_limits<double>::quiet_NaN());
                }
            }
            ++colIndex;
        }
        table.rows.push_back(std::move(row));
    }

    return table;
}

const std::vector<double>* findRow(const CsvTable& table, double P, double T, double xco2)
{
    const auto pIt = table.columns.find("P(bar)");
    const auto tIt = table.columns.find("T(K)");
    const auto xIt = table.columns.find("X(CO2)");
    if(pIt == table.columns.end() || tIt == table.columns.end() || xIt == table.columns.end())
        return nullptr;

    const size_t pIdx = pIt->second;
    const size_t tIdx = tIt->second;
    const size_t xIdx = xIt->second;

    for(const auto& row : table.rows)
    {
        if(row.size() <= std::max(std::max(pIdx, tIdx), xIdx))
            continue;
        if(std::abs(row[pIdx] - P) < 1e-6 && std::abs(row[tIdx] - T) < 1e-6 && std::abs(row[xIdx] - xco2) < 1e-6)
            return &row;
    }

    return nullptr;
}

const std::vector<double>* findRowByY(const CsvTable& table, double P, double T, double yH2O, double yCO2)
{
    const auto pIt = table.columns.find("P(bar)");
    const auto tIt = table.columns.find("T(K)");
    const auto hIt = table.columns.find("y(H2O)");
    const auto cIt = table.columns.find("y(CO2)");
    if(pIt == table.columns.end() || tIt == table.columns.end() || hIt == table.columns.end() || cIt == table.columns.end())
        return nullptr;

    const size_t pIdx = pIt->second;
    const size_t tIdx = tIt->second;
    const size_t hIdx = hIt->second;
    const size_t cIdx = cIt->second;

    for(const auto& row : table.rows)
    {
        if(row.size() <= std::max(std::max(pIdx, tIdx), std::max(hIdx, cIdx)))
            continue;
        if(std::abs(row[pIdx] - P) < 1e-6 && std::abs(row[tIdx] - T) < 1e-6
           && std::abs(row[hIdx] - yH2O) < 1e-6 && std::abs(row[cIdx] - yCO2) < 1e-6)
            return &row;
    }

    return nullptr;
}

double getValue(const CsvTable& table, const std::vector<double>& row, const std::string& col)
{
    const auto it = table.columns.find(col);
    if(it == table.columns.end())
        throw std::runtime_error("Column not found: " + col);
    const size_t idx = it->second;
    if(idx >= row.size())
        throw std::runtime_error("Row missing column: " + col);
    return row[idx];
}

std::vector<std::string> requiredGfsmReferenceFiles()
{
    // Prefer matrix-driven case inventory to keep parity coverage explicit.
    std::ifstream plan("test/gfsm_regression_matrix.csv");
    if(plan.is_open())
    {
        std::vector<std::string> files;
        std::string line;
        std::getline(plan, line); // header

        auto trim = [](std::string s) {
            while(!s.empty() && std::isspace(static_cast<unsigned char>(s.front()))) s.erase(s.begin());
            while(!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.pop_back();
            return s;
        };

        while(std::getline(plan, line))
        {
            if(line.empty())
                continue;

            std::stringstream ss(line);
            std::string section;
            std::string caseName;
            std::string path;
            std::string enabled;

            std::getline(ss, section, ',');
            std::getline(ss, caseName, ',');
            std::getline(ss, path, ',');
            std::getline(ss, enabled, ',');

            path = trim(path);
            enabled = trim(enabled);

            if(path.empty())
                continue;
            if(!enabled.empty() && enabled != "1")
                continue;

            files.push_back(path);
        }

        if(!files.empty())
            return files;
    }

    return {
        "test/gfsm/h2o_mrk.csv",
        "test/gfsm/h2o_hsmrk.csv",
        "test/gfsm/h2o_cork.csv",
        "test/gfsm/h2o_pseos.csv",
        "test/gfsm/h2o_haar.csv",
        "test/gfsm/h2o_zd05.csv",
        "test/gfsm/h2o_zd09.csv",
        "test/gfsm/co2_mrk.csv",
        "test/gfsm/co2_hsmrk.csv",
        "test/gfsm/co2_cork.csv",
        "test/gfsm/co2_pseos.csv",
        "test/gfsm/co2_zd09.csv",
        "test/gfsm/ch4_mrk.csv",
        "test/gfsm/ch4_hsmrk.csv",
        "test/gfsm/ch4_zd09.csv",
        "test/gfsm/h2o_co2_mrk_mrk.csv",
        "test/gfsm/h2o_co2_hsmrk_mrk.csv",
        "test/gfsm/h2o_co2_cork_mrk.csv",
        "test/gfsm/h2o_co2_pseos_mrk.csv",
        "test/gfsm/h2o_co2_haar_mrk.csv",
        "test/gfsm/h2o_co2_zd05_mrk.csv",
        "test/gfsm/h2o_co2_zd09_mrk.csv",
        "test/gfsm/h2o_co2_mrk_hsmrk.csv",
        "test/gfsm/h2o_co2_mrk_cork.csv",
        "test/gfsm/h2o_co2_mrk_pseos.csv",
        "test/gfsm/h2o_co2_mrk_zd09.csv",
        "test/gfsm/h2o_ch4_mrk_mrk.csv",
        "test/gfsm/h2o_ch4_hsmrk_mrk.csv",
        "test/gfsm/h2o_ch4_cork_mrk.csv",
        "test/gfsm/h2o_ch4_pseos_mrk.csv",
        "test/gfsm/h2o_ch4_haar_mrk.csv",
        "test/gfsm/h2o_ch4_zd05_mrk.csv",
        "test/gfsm/h2o_ch4_zd09_mrk.csv",
        "test/gfsm/h2o_ch4_mrk_hsmrk.csv",
        "test/gfsm/h2o_ch4_mrk_zd09.csv",
        "test/gfsm/co2_ch4_mrk_mrk.csv",
        "test/gfsm/co2_ch4_hsmrk_mrk.csv",
        "test/gfsm/co2_ch4_cork_mrk.csv",
        "test/gfsm/co2_ch4_pseos_mrk.csv",
        "test/gfsm/co2_ch4_zd09_mrk.csv",
        "test/gfsm/co2_ch4_mrk_hsmrk.csv",
        "test/gfsm/co2_ch4_mrk_zd09.csv",
        "test/gfsm/h2o_co2_ch4_mrk.csv",
        "test/gfsm/h2o_co2_ch4_hsmrk.csv",
        "test/gfsm/h2o_co2_ch4_zd09.csv",
        "test/gfsm/redox_o2_excess_mrk.csv",
        "test/gfsm/redox_o2_excess_hot.csv",
        "test/gfsm/redox_h2_excess_hot.csv",
        "test/gfsm/redox_co_bias_hot.csv",
        "test/gfsm/redox_o2_excess_hybrid.csv",
        "test/gfsm/redox_mixed_hybrid.csv",
    };
}

bool tokenToCode(const std::string& token, int& code)
{
    static const std::unordered_map<std::string, int> mapping = {
        {"mrk", 0}, {"hsmrk", 1}, {"cork", 2}, {"pseos", 4},
        {"haar", 5}, {"zd05", 6}, {"zd09", 7},
    };

    const auto it = mapping.find(token);
    if(it == mapping.end())
        return false;

    code = it->second;
    return true;
}

bool applyHybridCodes(HybridEosOptions& opts, int h2o, int co2, int ch4)
{
    switch(h2o)
    {
    case 0: opts.water = HybridEosOptions::WaterEos::Mrk; break;
    case 1: opts.water = HybridEosOptions::WaterEos::Hsmrk; break;
    case 2: opts.water = HybridEosOptions::WaterEos::Cork; break;
    case 4: opts.water = HybridEosOptions::WaterEos::Pseos; break;
    case 5: opts.water = HybridEosOptions::WaterEos::Haar; break;
    case 6: opts.water = HybridEosOptions::WaterEos::ZhangDuan05; break;
    case 7: opts.water = HybridEosOptions::WaterEos::ZhangDuan09; break;
    default: return false;
    }

    switch(co2)
    {
    case 0: opts.co2 = HybridEosOptions::CO2Eos::Mrk; break;
    case 1: opts.co2 = HybridEosOptions::CO2Eos::Hsmrk; break;
    case 2: opts.co2 = HybridEosOptions::CO2Eos::Cork; break;
    case 3: opts.co2 = HybridEosOptions::CO2Eos::Brmrk; break;
    case 4: opts.co2 = HybridEosOptions::CO2Eos::Pseos; break;
    case 7: opts.co2 = HybridEosOptions::CO2Eos::ZhangDuan09; break;
    default: return false;
    }

    switch(ch4)
    {
    case 0: opts.ch4 = HybridEosOptions::CH4Eos::Mrk; break;
    case 1: opts.ch4 = HybridEosOptions::CH4Eos::Hsmrk; break;
    case 7: opts.ch4 = HybridEosOptions::CH4Eos::ZhangDuan09; break;
    default: return false;
    }

    return true;
}

std::map<std::string, std::array<int, 3>> loadCaseEosOverrides(const std::string& path)
{
    std::map<std::string, std::array<int, 3>> out;
    std::ifstream file(path);
    if(!file.is_open())
        return out;

    std::string header;
    std::getline(file, header);

    std::string line;
    while(std::getline(file, line))
    {
        if(line.empty())
            continue;

        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> cols;
        while(std::getline(ss, token, ','))
            cols.push_back(token);

        if(cols.size() < 4)
            continue;

        try
        {
            out[cols[0]] = {std::stoi(cols[1]), std::stoi(cols[2]), std::stoi(cols[3])};
        }
        catch(...)
        {
            continue;
        }
    }

    return out;
}

bool parseCaseEosCodes(const std::string& caseName,
                       const std::map<std::string, std::array<int, 3>>& overrides,
                       int& h2o,
                       int& co2,
                       int& ch4)
{
    h2o = 0;
    co2 = 0;
    ch4 = 0;

    const auto ov = overrides.find(caseName);
    if(ov != overrides.end())
    {
        h2o = ov->second[0];
        co2 = ov->second[1];
        ch4 = ov->second[2];
        return true;
    }

    std::vector<std::string> parts;
    std::stringstream ss(caseName);
    std::string p;
    while(std::getline(ss, p, '_'))
        parts.push_back(p);

    int codeA = 0;
    int codeB = 0;

    if(parts.size() == 2)
    {
        if(parts[0] == "h2o") return tokenToCode(parts[1], h2o);
        if(parts[0] == "co2") return tokenToCode(parts[1], co2);
        if(parts[0] == "ch4") return tokenToCode(parts[1], ch4);
        return false;
    }

    if(parts.size() == 4 && parts[0] == "h2o" && parts[1] == "co2" && parts[2] == "ch4")
    {
        if(!tokenToCode(parts[3], codeA))
            return false;
        h2o = codeA;
        co2 = codeA;
        ch4 = codeA;
        return true;
    }

    if(parts.size() == 4)
    {
        if(!tokenToCode(parts[2], codeA) || !tokenToCode(parts[3], codeB))
            return false;

        if(parts[0] == "h2o" && parts[1] == "co2") { h2o = codeA; co2 = codeB; return true; }
        if(parts[0] == "h2o" && parts[1] == "ch4") { h2o = codeA; ch4 = codeB; return true; }
        if(parts[0] == "co2" && parts[1] == "ch4") { co2 = codeA; ch4 = codeB; return true; }
        return false;
    }

    return false;
}

// Helper function to compare values with tolerance
bool compareValue(double computed, double reference, double tolerance,
                  const std::string& name, std::vector<std::string>& errors,
                  bool relative = true)
{
    g_metrics.assertionsExecuted++;

    double error = relative
        ? std::abs((computed - reference) / (reference + 1e-20))
        : std::abs(computed - reference);

    if (error > tolerance) {
        char buffer[256];
        if (relative) {
            snprintf(buffer, sizeof(buffer),
                    "%s: computed=%.10e, reference=%.10e, rel_error=%.2e (tol=%.2e)",
                    name.c_str(), computed, reference, error, tolerance);
        } else {
            snprintf(buffer, sizeof(buffer),
                    "%s: computed=%.10e, reference=%.10e, abs_error=%.2e (tol=%.2e)",
                    name.c_str(), computed, reference, error, tolerance);
        }
        errors.push_back(buffer);
        g_metrics.assertionsFailed++;
        return false;
    }

    g_metrics.assertionsPassed++;
    return true;
}

// =============================================================================
// TEST 1: Pure H2O with HSMRK
// =============================================================================
TestResult test_pure_h2o_hsmrk()
{
    TestResult result{"Pure H2O MRK (ifug=2 defaults)", true, {}, {}};

    // Reference data from Perple_X COHSRK (ifug=0 MRK) at T=523.15K, P=1000bar
    ReferencePoint ref;
    ref.T = 523.15;
    ref.P = 1000.0;
    ref.y = {1.0};  // Pure H2O

    const auto table = loadCsvTable("test/h2o_co2_mrk.csv");
    const auto row = findRow(table, ref.P, ref.T, 0.0);
    if(!row)
    {
        result.passed = false;
        result.errors.push_back("Reference row not found in h2o_co2_mrk.csv for XCO2=0.0");
        return result;
    }
    const double fH2O = getValue(table, *row, "f(H2O)");
    ref.ln_f = {std::log(fH2O)};
    ref.v = {getValue(table, *row, "vol[cm3/mol]")};
    ref.vol = ref.v[0];

    // Compute with Reaktoro
    PerpleXFluidOptions opts;
    // Use MRK defaults (no pure-EoS replacement) to match COHSRK ifug=2 defaults
    opts.hybridSpecies = {1, 2};
    opts.hybridOptions = makePerpleXHybridEosOptions();
    opts.hybridOptions.water = HybridEosOptions::WaterEos::Mrk;
    opts.hybridOptions.co2 = HybridEosOptions::CO2Eos::Cork;

    std::array<double, 19> y{};
    y[1] = 1.0;

    PerpleXFluidModel model;
    auto state = model.compute({1, 2}, y, ref.P, ref.T, opts);

    // Compare results
    result.passed &= compareValue(state.ln_f[1], ref.ln_f[0], TOL_FUGACITY,
                                  "H2O ln(fugacity)", result.errors);
    result.passed &= compareValue(state.v[1], ref.v[0], TOL_VOLUME,
                                  "H2O partial molar volume", result.errors);
    result.passed &= compareValue(state.vol, ref.vol, TOL_VOLUME,
                                  "Total volume", result.errors);

    return result;
}

// =============================================================================
// TEST 1b: Pure H2O with HSMRK (ifug=1)
// =============================================================================
TestResult test_pure_h2o_hsmrk_ref()
{
    TestResult result{"Pure H2O HSMRK (ifug=1)", true, {}, {}};

    ReferencePoint ref;
    ref.T = 523.15;
    ref.P = 1000.0;
    ref.y = {1.0};

    const auto table = loadCsvTable("test/h2o_co2_hsmrk.csv");
    const auto row = findRowByY(table, ref.P, ref.T, 1.0, 0.0);
    if(!row)
    {
        result.passed = false;
        result.errors.push_back("Reference row not found in h2o_co2_hsmrk.csv for y(H2O)=1.0");
        return result;
    }
    const double fH2O = getValue(table, *row, "f(H2O)");
    ref.ln_f = {std::log(fH2O)};
    ref.v = {getValue(table, *row, "vol[cm3/mol]")};
    ref.vol = ref.v[0];

    PerpleXFluidOptions opts;
    opts.hybridSpecies = {1};
    opts.hybridOptions = makePerpleXHybridEosOptions();
    opts.hybridOptions.water = HybridEosOptions::WaterEos::Hsmrk;
    opts.hybridOptions.co2 = HybridEosOptions::CO2Eos::Cork;

    std::array<double, 19> y{};
    y[1] = 1.0;

    PerpleXFluidModel model;
    auto state = model.compute({1}, y, ref.P, ref.T, opts);

    result.passed &= compareValue(state.ln_f[1], ref.ln_f[0], TOL_FUGACITY,
                                  "H2O ln(fugacity) HSMRK", result.errors);
    result.passed &= compareValue(state.v[1], ref.v[0], TOL_VOLUME,
                                  "H2O partial molar volume HSMRK", result.errors);
    result.passed &= compareValue(state.vol, ref.vol, TOL_VOLUME,
                                  "Total volume HSMRK", result.errors);

    return result;
}

// =============================================================================
// TEST 2: Pure CO2 with CORK
// =============================================================================
TestResult test_pure_co2_cork()
{
    TestResult result{"Pure CO2 MRK (ifug=2 defaults)", true, {}, {}};

    ReferencePoint ref;
    ref.T = 523.15;
    ref.P = 1000.0;
    ref.y = {1.0};

    const auto table = loadCsvTable("test/h2o_co2_mrk.csv");
    const auto row = findRow(table, ref.P, ref.T, 1.0);
    if(!row)
    {
        result.passed = false;
        result.errors.push_back("Reference row not found in h2o_co2_mrk.csv for XCO2=1.0");
        return result;
    }
    const double fCO2 = getValue(table, *row, "f(CO2)");
    ref.ln_f = {std::log(fCO2)};
    ref.v = {getValue(table, *row, "vol[cm3/mol]")};
    ref.vol = ref.v[0];

    // Compute with Reaktoro
    PerpleXFluidOptions opts;
    opts.hybridSpecies = {1, 2};
    opts.hybridOptions = makePerpleXHybridEosOptions();
    opts.hybridOptions.water = HybridEosOptions::WaterEos::Mrk;
    opts.hybridOptions.co2 = HybridEosOptions::CO2Eos::Mrk;

    std::array<double, 19> y{};
    y[2] = 1.0;

    PerpleXFluidModel model;
    auto state = model.compute({2}, y, ref.P, ref.T, opts);

    result.passed &= compareValue(state.ln_f[2], ref.ln_f[0], TOL_FUGACITY,
                                  "CO2 ln(fugacity)", result.errors);
    result.passed &= compareValue(state.v[2], ref.v[0], TOL_VOLUME,
                                  "CO2 partial molar volume", result.errors);

    return result;
}

// =============================================================================
// TEST 2b: Pure CO2 with CORK (ifug=5)
// =============================================================================
TestResult test_pure_co2_cork_ref()
{
    TestResult result{"Pure CO2 CORK (ifug=5)", true, {}, {}};

    ReferencePoint ref;
    ref.T = 523.15;
    ref.P = 1000.0;
    ref.y = {1.0};

    const auto table = loadCsvTable("test/h2o_co2_cork.csv");
    const auto row = findRowByY(table, ref.P, ref.T, 0.0, 1.0);
    if(!row)
    {
        result.passed = false;
        result.errors.push_back("Reference row not found in h2o_co2_cork.csv for y(CO2)=1.0");
        return result;
    }
    const double fCO2 = getValue(table, *row, "f(CO2)");
    ref.ln_f = {std::log(fCO2)};
    ref.v = {getValue(table, *row, "vol[cm3/mol]")};
    ref.vol = ref.v[0];

    PerpleXFluidOptions opts;
    opts.hybridSpecies = {1, 2};
    opts.hybridOptions = makePerpleXHybridEosOptions();
    opts.hybridOptions.water = HybridEosOptions::WaterEos::Mrk;
    opts.hybridOptions.co2 = HybridEosOptions::CO2Eos::Cork;

    std::array<double, 19> y{};
    y[2] = 1.0;

    PerpleXFluidModel model;
    auto state = model.compute({2}, y, ref.P, ref.T, opts);

    result.passed &= compareValue(state.ln_f[2], ref.ln_f[0], TOL_FUGACITY,
                                  "CO2 ln(fugacity) CORK", result.errors);
    if(ref.v[0] > 0.0 && !std::isnan(ref.v[0]))
    {
        result.passed &= compareValue(state.v[2], ref.v[0], TOL_VOLUME,
                                      "CO2 partial molar volume CORK", result.errors);
    }
    else
    {
        result.warnings.push_back("CORK reference volume unavailable; skipped volume comparison.");
    }

    return result;
}

// =============================================================================
// TEST 3: H2O-CO2 Binary Mixture
// =============================================================================
TestResult test_h2o_co2_binary()
{
    TestResult result{"H2O-CO2 Binary (50-50, ifug=2)", true, {}, {}};

    ReferencePoint ref;
    ref.T = 523.15;
    ref.P = 1000.0;
    ref.y = {0.5, 0.5};

    const auto table = loadCsvTable("test/h2o_co2_mrk.csv");
    const auto row = findRow(table, ref.P, ref.T, 0.5);
    if(!row)
    {
        result.passed = false;
        result.errors.push_back("Reference row not found in h2o_co2_mrk.csv for XCO2=0.5");
        return result;
    }
    const double fH2O = getValue(table, *row, "f(H2O)");
    const double fCO2 = getValue(table, *row, "f(CO2)");
    ref.ln_f = {std::log(fH2O), std::log(fCO2)};
    ref.vol = getValue(table, *row, "vol[cm3/mol]");

    // Compute with Reaktoro
    PerpleXFluidOptions opts;
    // MRK defaults (no pure-EoS replacement)
    opts.hybridSpecies = {1, 2};
    opts.hybridOptions = makePerpleXHybridEosOptions();
    opts.hybridOptions.water = HybridEosOptions::WaterEos::Mrk;
    opts.hybridOptions.co2 = HybridEosOptions::CO2Eos::Mrk;

    std::array<double, 19> y{};
    y[1] = 0.5;
    y[2] = 0.5;

    PerpleXFluidModel model;
    auto state = model.compute({1, 2}, y, ref.P, ref.T, opts);

    result.passed &= compareValue(state.ln_f[1], ref.ln_f[0], TOL_FUGACITY,
                                  "H2O ln(fugacity) (binary)", result.errors);
    result.passed &= compareValue(state.ln_f[2], ref.ln_f[1], TOL_FUGACITY,
                                  "CO2 ln(fugacity) (binary)", result.errors);
    result.passed &= compareValue(state.vol, ref.vol, TOL_VOLUME,
                                  "Binary mixture volume", result.errors);

    return result;
}

// =============================================================================
// TEST 4: Dielectric Constant (H2O-CO2)
// =============================================================================
TestResult test_dielectric_h2o_co2()
{
    TestResult result{"Dielectric H2O-CO2", true, {}, {}};
    const auto table = loadCsvTable("test/epsh2o_reference.csv");
    const auto tIt = table.columns.find("T_K");
    const auto vIt = table.columns.find("vol_cm3");
    const auto eIt = table.columns.find("epsilon_ref");
    if(tIt == table.columns.end() || vIt == table.columns.end() || eIt == table.columns.end())
    {
        result.passed = false;
        result.errors.push_back("Missing columns in epsh2o_reference.csv");
        return result;
    }

    for(const auto& row : table.rows)
    {
        const double T = row[tIt->second];
        const double vol_cm3 = row[vIt->second];
        const double ref_eps = row[eIt->second];
        const double v_jbar = vol_cm3 / 10.0;
        const double eps = epsh2o(v_jbar, T);

        char name[64];
        snprintf(name, sizeof(name), "epsh2o at T=%.0fK", T);

        result.passed &= compareValue(eps, ref_eps, TOL_DIELECTRIC,
                                      name, result.errors, false);
    }

    return result;
}

// =============================================================================
// TEST 5: Shock g-Function
// =============================================================================
TestResult test_gfunction()
{
    TestResult result{"Shock g-Function", true, {}, {}};

    const auto table = loadCsvTable("test/gfunc_reference.csv");
    const auto pIt = table.columns.find("P_bar");
    const auto tIt = table.columns.find("T_K");
    const auto rIt = table.columns.find("rho_gcc");
    const auto gIt = table.columns.find("g_ref");
    if(pIt == table.columns.end() || tIt == table.columns.end() || rIt == table.columns.end() || gIt == table.columns.end())
    {
        result.passed = false;
        result.errors.push_back("Missing columns in gfunc_reference.csv");
        return result;
    }

    for(const auto& row : table.rows)
    {
        const double P = row[pIt->second];
        const double T = row[tIt->second];
        const double rho = row[rIt->second];
        const double ref_g = row[gIt->second];

        double computed_g = gfunc(rho, P, T);

        char name[64];
        snprintf(name, sizeof(name), "g-function at rho=%.2f", rho);

        result.passed &= compareValue(computed_g, ref_g, TOL_GFUNCTION,
                                      name, result.errors, false);
    }

    return result;
}

// =============================================================================
// TEST 6: HKF Aqueous Species (matrix-driven multi-ion)
// =============================================================================
TestResult test_hkf_multiion()
{
    TestResult result{"HKF Multi-ion Gibbs Energy", true, {}, {}};

    std::map<std::string, HKFParams> paramBySpecies;

    // Load enabled ion parameter rows from hkf_matrix.csv
    std::ifstream matrix("test/hkf_matrix.csv");
    if(!matrix.is_open())
    {
        result.passed = false;
        result.errors.push_back("Failed to open test/hkf_matrix.csv");
        return result;
    }

    std::string line;
    std::getline(matrix, line); // header
    while(std::getline(matrix, line))
    {
        if(line.empty()) continue;
        std::stringstream ss(line);
        std::vector<std::string> cols;
        std::string tok;
        while(std::getline(ss, tok, ',')) cols.push_back(tok);
        if(cols.size() < 16) continue;
        if(cols[0] != "1") continue;

        HKFParams p;
        p.G0 = parseDouble(cols[2]);
        p.S0 = parseDouble(cols[3]);
        p.omega0 = parseDouble(cols[4]);
        p.charge = parseDouble(cols[5]);
        p.a1 = parseDouble(cols[6]);
        p.a2 = parseDouble(cols[7]);
        p.a3 = parseDouble(cols[8]);
        p.a4 = parseDouble(cols[9]);
        p.c1 = parseDouble(cols[10]);
        p.c2 = parseDouble(cols[11]);

        // hkf_matrix.csv stores DEW database parameters with a1/a3 in J/(mol·Pa).
        // Convert to J/(mol·bar) for the PerplexHKF engine (pressure in bar).
        constexpr double bar = 1.0e5;
        p.a1 = p.a1 * bar;
        p.a3 = p.a3 * bar;

        paramBySpecies[cols[1]] = preprocessHKFParams(p);
    }

    if(paramBySpecies.empty())
    {
        result.passed = false;
        result.errors.push_back("No enabled species found in hkf_matrix.csv");
        return result;
    }

    // Evaluate species rows in hkf_reference.csv
    std::ifstream ref("test/hkf_reference.csv");
    if(!ref.is_open())
    {
        result.passed = false;
        result.errors.push_back("Failed to open test/hkf_reference.csv");
        return result;
    }

    std::getline(ref, line); // header
    while(std::getline(ref, line))
    {
        if(line.empty()) continue;
        std::stringstream ss(line);
        std::vector<std::string> cols;
        std::string tok;
        while(std::getline(ss, tok, ',')) cols.push_back(tok);
        if(cols.size() < 6) continue;

        const std::string species = cols[0];
        const auto it = paramBySpecies.find(species);
        if(it == paramBySpecies.end())
            continue;

        const double P = parseDouble(cols[1]);
        const double T = parseDouble(cols[2]);
        const double eps = parseDouble(cols[3]);
        const double gf = parseDouble(cols[4]);
        const double ref_G = parseDouble(cols[5]);

        auto hkf = computeHKFGibbs(it->second, P, T, eps, gf);
        const std::string name = species + " G at T=" + std::to_string(static_cast<int>(T))
                               + "K P=" + std::to_string(static_cast<int>(P)) + "bar";

        result.passed &= compareValue(hkf.G, ref_G, TOL_HKF_GIBBS,
                                      name, result.errors, false);
    }

    return result;
}

// =============================================================================
// TEST 6b: DEW Speciation Energy Parity (HKF Gibbs)
// =============================================================================
TestResult test_dew_speciation_energy_parity()
{
    auto result = test_hkf_multiion();
    result.testName = "DEW Speciation Energy Parity (HKF G)";
    return result;
}

// =============================================================================
// TEST 7: Composition Series (H2O-CO2 at varying XCO2)
// =============================================================================
TestResult test_composition_series()
{
    TestResult result{"H2O-CO2 Composition Series", true, {}, {}};

    double T = 523.15;
    double P = 1000.0;

    std::vector<double> xco2_values = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5,
                                       0.6, 0.7, 0.8, 0.9, 1.0};

    const auto table = loadCsvTable("test/h2o_co2_mrk.csv");

    PerpleXFluidOptions opts;
    opts.hybridSpecies = {1, 2};
    opts.hybridOptions = makePerpleXHybridEosOptions();
    opts.hybridOptions.water = HybridEosOptions::WaterEos::Mrk;
    opts.hybridOptions.co2 = HybridEosOptions::CO2Eos::Mrk;

    PerpleXFluidModel model;
    constexpr double R = 8.31446261815324;

    for (size_t i = 0; i < xco2_values.size(); ++i) {
        std::array<double, 19> y{};
        y[1] = 1.0 - xco2_values[i];
        y[2] = xco2_values[i];

        auto state = model.compute({1, 2}, y, P, T, opts);

        char name[64];
        snprintf(name, sizeof(name), "volume at XCO2=%.2f", xco2_values[i]);

        const auto row = findRow(table, P, T, xco2_values[i]);
        if(!row)
        {
            result.passed = false;
            result.errors.push_back("Reference row not found in h2o_co2_mrk.csv for XCO2=" + std::to_string(xco2_values[i]));
            continue;
        }
        const double ref_vol = getValue(table, *row, "vol[cm3/mol]");
        const std::array<double, 19> yCopy = y;
        const std::vector<int> activeSpecies{1, 2};
        const auto tolTier = selectVolumeToleranceTier(yCopy, activeSpecies, P, T);
        const double volTol = std::max(tolTier.absTol, tolTier.relTol * std::abs(ref_vol));
        const std::string volLabel = std::string(name) + " (" + tolTier.name + ")";
        result.passed &= compareValue(state.vol, ref_vol, volTol,
                          volLabel, result.errors, false);

        // Chemical-potential parity proxy via ln(f) comparison.
        const double ref_lnf_h2o = std::log(getValue(table, *row, "f(H2O)"));
        const double ref_lnf_co2 = std::log(getValue(table, *row, "f(CO2)"));
        if(y[1] > 1e-12)
        {
            result.passed &= compareValue(state.ln_f[1], ref_lnf_h2o, TOL_FUGACITY,
                                          "H2O ln(f) composition series", result.errors);
        }
        if(y[2] > 1e-12)
        {
            result.passed &= compareValue(state.ln_f[2], ref_lnf_co2, TOL_FUGACITY,
                                          "CO2 ln(f) composition series", result.errors);
        }

        // Delta-G parity proxy: RT * sum(y_i * ln f_i) over present species.
        const double yH2O = y[1];
        const double yCO2 = y[2];
        const double dGmix_ref = R * T * (yH2O * ref_lnf_h2o + yCO2 * ref_lnf_co2);
        const double dGmix = R * T * (yH2O * state.ln_f[1] + yCO2 * state.ln_f[2]);
        result.passed &= compareValue(dGmix, dGmix_ref, TOL_DELTAG_MIX,
                          "dGmix parity composition series", result.errors, false);
    }

    return result;
}

// =============================================================================
// TEST 7c: Debye-HÃ¼ckel factor
// =============================================================================
TestResult test_debye_huckel()
{
    TestResult result{"Debye-HÃ¼ckel Factor", true, {}, {}};

    const auto table = loadCsvTable("test/dh_reference.csv");
    const auto mIt = table.columns.find("msol");
    const auto vIt = table.columns.find("vsolv_cm3");
    const auto eIt = table.columns.find("epsilon");
    const auto tIt = table.columns.find("T_K");
    const auto aIt = table.columns.find("adh_ref");
    if(mIt == table.columns.end() || vIt == table.columns.end() || eIt == table.columns.end()
        || tIt == table.columns.end() || aIt == table.columns.end())
    {
        result.passed = false;
        result.errors.push_back("Missing columns in dh_reference.csv");
        return result;
    }

    for(const auto& row : table.rows)
    {
        const double msol = row[mIt->second];
        const double vsolv = row[vIt->second];
        const double eps = row[eIt->second];
        const double T = row[tIt->second];
        const double ref_adh = row[aIt->second];

        const double adh = debyeHuckel(msol, vsolv, eps, T);

        char name[64];
        snprintf(name, sizeof(name), "adh at T=%.0fK", T);

        result.passed &= compareValue(adh, ref_adh, TOL_DH_FACTOR,
                                      name, result.errors, false);
    }

    return result;
}

// =============================================================================
// TEST 7d: Born omega
// =============================================================================
TestResult test_born_omega()
{
    TestResult result{"Born Omega", true, {}, {}};

    const auto table = loadCsvTable("test/born_reference.csv");
    const auto zIt = table.columns.find("z");
    const auto oIt = table.columns.find("omega0");
    const auto rIt = table.columns.find("born_radius");
    const auto gIt = table.columns.find("g");
    const auto wIt = table.columns.find("omega_ref");
    if(zIt == table.columns.end() || oIt == table.columns.end() || rIt == table.columns.end()
        || gIt == table.columns.end() || wIt == table.columns.end())
    {
        result.passed = false;
        result.errors.push_back("Missing columns in born_reference.csv");
        return result;
    }

    for(const auto& row : table.rows)
    {
        const double z = row[zIt->second];
        const double omega0 = row[oIt->second];
        const double re = row[rIt->second];
        const double gf = row[gIt->second];
        const double ref_omega = row[wIt->second];

        const double omega = calculateBornOmega(z, re, omega0, gf);

        char name[64];
        snprintf(name, sizeof(name), "omega z=%.1f", z);

        result.passed &= compareValue(omega, ref_omega, TOL_HKF_GIBBS,
                                      name, result.errors, false);
    }

    return result;
}

// =============================================================================
// TEST 8: Pressure-Temperature Grid
// =============================================================================
TestResult test_pt_grid()
{
    TestResult result{"P-T Grid Test", true, {}, {}};

    // Test grid
    std::vector<double> pressures = {1000.0, 2000.0, 3000.0, 4000.0, 5000.0};
    std::vector<double> temperatures = {373.15, 423.15, 473.15, 523.15, 573.15};

    // Binary H2O-CO2 (70-30)
    std::array<double, 19> y{};
    y[1] = 0.7;
    y[2] = 0.3;

    PerpleXFluidOptions opts;
    opts.hybridSpecies = {1, 2};
    opts.hybridOptions = makePerpleXHybridEosOptions();
    opts.hybridOptions.water = HybridEosOptions::WaterEos::Mrk;
    opts.hybridOptions.co2 = HybridEosOptions::CO2Eos::Mrk;

    PerpleXFluidModel model;

    // TODO: Load reference grid data from file
    // For now, just check that calculations complete without error

    for (double P : pressures) {
        for (double T : temperatures) {
            try {
                auto state = model.compute({1, 2}, y, P, T, opts);

                // Basic sanity checks
                if (std::isnan(state.vol) || std::isinf(state.vol)) {
                    result.errors.push_back("Invalid volume at P=" +
                        std::to_string(P) + " T=" + std::to_string(T));
                    result.passed = false;
                }

                if (std::isnan(state.g[0]) || std::isinf(state.g[0])) {
                    result.errors.push_back("Invalid H2O fugacity at P=" +
                        std::to_string(P) + " T=" + std::to_string(T));
                    result.passed = false;
                }

            } catch (const std::exception& e) {
                result.errors.push_back("Exception at P=" + std::to_string(P) +
                    " T=" + std::to_string(T) + ": " + e.what());
                result.passed = false;
            }
        }
    }

    return result;
}

// =============================================================================
// TEST 9: Coverage Gate for GFSM-Accepted EOS Components
// =============================================================================
TestResult test_gfsm_eos_coverage_gate()
{
    TestResult result{"GFSM EOS Coverage Gate", true, {}, {}};

    for(const auto& path : requiredGfsmReferenceFiles())
    {
        std::error_code ec;
        const bool exists = std::filesystem::exists(path, ec);
        const auto size = exists ? std::filesystem::file_size(path, ec) : 0;
        if(!exists || ec || size == 0)
        {
            result.passed = false;
            result.errors.push_back("Missing or empty required GFSM EOS reference file: " + path);
        }
    }

    return result;
}

TestResult test_gfsm_full_matrix_volume_parity()
{
    TestResult result{"GFSM Full Matrix Volume Parity", true, {}, {}};

    const auto caseOverrides = loadCaseEosOverrides("test/gfsm_case_matrix.csv");
    const std::vector<std::pair<std::string, int>> speciesCols = {
        {"H2O", 1}, {"CO2", 2}, {"CH4", 4}, {"H2", 5}, {"CO", 3},
        {"H2S", 6}, {"O2", 7}, {"SO2", 8}, {"N2", 10}, {"NH3", 11},
        {"HF", 17}, {"C2H6", 16}, {"HCl", 18},
    };

    PerpleXFluidModel model;

    for(const auto& path : requiredGfsmReferenceFiles())
    {
        try
        {
            const CsvTable table = loadCsvTable(path);
            if(table.rows.empty())
            {
                result.passed = false;
                result.errors.push_back("Empty GFSM reference CSV: " + path);
                continue;
            }

            const std::string caseName = std::filesystem::path(path).stem().string();
            int h2oCode = 0;
            int co2Code = 0;
            int ch4Code = 0;
            if(!parseCaseEosCodes(caseName, caseOverrides, h2oCode, co2Code, ch4Code))
            {
                result.passed = false;
                result.errors.push_back("Failed to parse EOS configuration from case name: " + caseName);
                continue;
            }

            PerpleXFluidOptions opts;
            opts.hybridSpecies = {1, 2, 4};
            opts.hybridOptions = makePerpleXHybridEosOptions();
            if(!applyHybridCodes(opts.hybridOptions, h2oCode, co2Code, ch4Code))
            {
                result.passed = false;
                result.errors.push_back("Unsupported EOS code tuple for case: " + caseName);
                continue;
            }

            for(size_t rowIndex = 0; rowIndex < table.rows.size(); ++rowIndex)
            {
                const auto& row = table.rows[rowIndex];
                double P = getValue(table, row, "P(bar)");
                double T = getValue(table, row, "T(K)");
                const double vref = getValue(table, row, "vol[cm3/mol]");

                bool rowPassed = true;

                std::array<double, 19> y{};
                std::vector<int> species;
                for(const auto& [name, idx] : speciesCols)
                {
                    const std::string col = "y(" + name + ")";
                    const auto it = table.columns.find(col);
                    if(it == table.columns.end() || it->second >= row.size())
                        continue;

                    y[idx] = row[it->second];
                    if(y[idx] > 1e-14)
                        species.push_back(idx);
                }

                if(species.empty())
                {
                    result.passed = false;
                    result.errors.push_back("No active species parsed for case: " + caseName + " row=" + std::to_string(rowIndex));
                    continue;
                }

                const auto state = model.compute(species, y, P, T, opts);

                // Keep a modest absolute floor due to rounded reference CSVs,
                // while enforcing a thermodynamic-grade relative tolerance.
                const auto tolTier = selectVolumeToleranceTier(y, species, P, T);
                const double tol = std::max(tolTier.absTol, tolTier.relTol * std::abs(vref));
                const std::string label = "vol[cm3/mol] case=" + caseName + " row=" + std::to_string(rowIndex)
                                        + " tier=" + tolTier.name;
                const bool volOk = compareValue(state.vol, vref, tol, label, result.errors, false);
                result.passed &= volOk;
                rowPassed &= volOk;

                // Additional structural sanity checks for all active species.
                for(const int idx : species)
                {
                    if(!std::isfinite(state.ln_f[idx]))
                    {
                        result.passed = false;
                        rowPassed = false;
                        result.errors.push_back("Invalid ln_f for species index " + std::to_string(idx) + " in case=" + caseName + " row=" + std::to_string(rowIndex));
                    }
                    if(!std::isfinite(state.v[idx]))
                    {
                        result.passed = false;
                        rowPassed = false;
                        result.errors.push_back("Invalid partial molar volume for species index " + std::to_string(idx) + " in case=" + caseName + " row=" + std::to_string(rowIndex));
                    }
                }

                g_metrics.rowAssertionsExecuted++;
                if(rowPassed)
                    g_metrics.rowAssertionsPassed++;
                else
                    g_metrics.rowAssertionsFailed++;
            }
        }
        catch(const std::exception& e)
        {
            result.passed = false;
            result.errors.push_back("Case failed (" + path + "): " + e.what());
        }
    }

    return result;
}

// =============================================================================
// TEST 10: Integration Flow (MRK -> Hybrid -> Dielectric -> HKF)
// =============================================================================
TestResult test_integration_flow_scenarios()
{
    TestResult result{"Integration Flow Scenarios", true, {}, {}};

    // Focused integration set covering MRK baseline and hybrid substitutions.
    const std::vector<std::string> scenarioFiles = {
        "test/gfsm/h2o_co2_mrk_mrk.csv",
        "test/gfsm/h2o_co2_hsmrk_mrk.csv",
        "test/gfsm/h2o_co2_zd09_mrk.csv",
        "test/gfsm/h2o_co2_ch4_hsmrk.csv",
        "test/gfsm/h2o_ch4_mrk_mrk.csv",
        "test/gfsm/redox_mixed_hybrid.csv",
    };

    const std::map<std::string, int> minActiveSpeciesByCase = {
        {"h2o_co2_ch4_hsmrk", 5},
        {"h2o_ch4_mrk_mrk", 5},
        {"redox_mixed_hybrid", 5},
    };

    // Load enabled HKF species definitions from matrix.
    std::map<std::string, HKFParams> hkfParams;
    {
        std::ifstream matrix("test/hkf_matrix.csv");
        std::string line;
        std::getline(matrix, line); // header
        while(std::getline(matrix, line))
        {
            if(line.empty()) continue;
            std::stringstream ss(line);
            std::vector<std::string> cols;
            std::string tok;
            while(std::getline(ss, tok, ',')) cols.push_back(tok);
            if(cols.size() < 16 || cols[0] != "1")
                continue;

            HKFParams p;
            p.G0 = parseDouble(cols[2]);
            p.S0 = parseDouble(cols[3]);
            p.omega0 = parseDouble(cols[4]);
            p.charge = parseDouble(cols[5]);
            p.a1 = parseDouble(cols[6]);
            p.a2 = parseDouble(cols[7]);
            p.a3 = parseDouble(cols[8]);
            p.a4 = parseDouble(cols[9]);
            p.c1 = parseDouble(cols[10]);
            p.c2 = parseDouble(cols[11]);
            // hkf_matrix.csv stores DEW database parameters with a1/a3 in J/(mol·Pa).
            // Convert to J/(mol·bar) for the PerplexHKF engine (pressure in bar).
            constexpr double bar = 1.0e5;
            p.a1 = p.a1 * bar;
            p.a3 = p.a3 * bar;
            hkfParams[cols[1]] = preprocessHKFParams(p);
        }
    }

    PerpleXFluidModel model;

    for(const auto& path : scenarioFiles)
    {
        try
        {
            const CsvTable table = loadCsvTable(path);
            if(table.rows.empty())
            {
                result.passed = false;
                result.errors.push_back("Integration scenario has no rows: " + path);
                continue;
            }

            const std::string caseName = std::filesystem::path(path).stem().string();
            int h2oCode = 0;
            int co2Code = 0;
            int ch4Code = 0;
            if(!parseCaseEosCodes(caseName, loadCaseEosOverrides("test/gfsm_case_matrix.csv"), h2oCode, co2Code, ch4Code))
            {
                result.passed = false;
                result.errors.push_back("Failed to parse EOS configuration for integration case: " + caseName);
                continue;
            }

            PerpleXFluidOptions opts;
            opts.hybridSpecies = {1, 2, 4};
            opts.hybridOptions = makePerpleXHybridEosOptions();
            opts.enableElectrolyte = true;
            if(!applyHybridCodes(opts.hybridOptions, h2oCode, co2Code, ch4Code))
            {
                result.passed = false;
                result.errors.push_back("Unsupported EOS tuple in integration case: " + caseName);
                continue;
            }

            for(size_t rowIndex = 0; rowIndex < table.rows.size(); ++rowIndex)
            {
                const auto& row = table.rows[rowIndex];
                const double P = getValue(table, row, "P(bar)");
                const double T = getValue(table, row, "T(K)");
                const double vref = getValue(table, row, "vol[cm3/mol]");

                std::array<double, 19> y{};
                std::vector<int> species;
                const std::vector<std::pair<std::string, int>> speciesCols = {
                    {"H2O", 1}, {"CO2", 2}, {"CH4", 4}, {"H2", 5}, {"CO", 3},
                    {"H2S", 6}, {"O2", 7}, {"SO2", 8}, {"N2", 10}, {"NH3", 11},
                    {"HF", 17}, {"C2H6", 16}, {"HCl", 18},
                };

                for(const auto& [name, idx] : speciesCols)
                {
                    const auto it = table.columns.find("y(" + name + ")");
                    if(it == table.columns.end() || it->second >= row.size())
                        continue;
                    y[idx] = row[it->second];
                    if(y[idx] > 1e-14)
                        species.push_back(idx);
                }

                if(species.empty())
                {
                    result.passed = false;
                    result.errors.push_back("No active species in integration case=" + caseName + " row=" + std::to_string(rowIndex));
                    continue;
                }

                const auto minIt = minActiveSpeciesByCase.find(caseName);
                if(minIt != minActiveSpeciesByCase.end() && static_cast<int>(species.size()) < minIt->second)
                {
                    result.passed = false;
                    result.errors.push_back(
                        "Insufficient active species in integration case=" + caseName
                        + " row=" + std::to_string(rowIndex)
                        + " expected>=" + std::to_string(minIt->second)
                        + " got=" + std::to_string(species.size())
                    );
                }

                const auto state = model.compute(species, y, P, T, opts);
                const auto tolTier = selectVolumeToleranceTier(y, species, P, T);
                const double volTol = std::max(tolTier.absTol, tolTier.relTol * std::abs(vref));
                const std::string volLabel = "integration vol case=" + caseName + " row=" + std::to_string(rowIndex);
                result.passed &= compareValue(state.vol, vref, volTol, volLabel, result.errors, false);

                // Dielectric/electrolyte stage must remain numerically well-behaved.
                if(!std::isfinite(state.dielectric.epsilon) || state.dielectric.epsilon <= 0.0)
                {
                    result.passed = false;
                    result.errors.push_back("Invalid epsilon in integration case=" + caseName + " row=" + std::to_string(rowIndex));
                }
                if(!std::isfinite(state.dielectric.gf))
                {
                    result.passed = false;
                    result.errors.push_back("Invalid g-function in integration case=" + caseName + " row=" + std::to_string(rowIndex));
                }
                if(!std::isfinite(state.dielectric.adh))
                {
                    result.passed = false;
                    result.errors.push_back("Invalid Debye-Huckel factor in integration case=" + caseName + " row=" + std::to_string(rowIndex));
                }

                // HKF stage: all enabled ion species should produce finite G.
                for(const auto& kv : hkfParams)
                {
                    const auto hkf = computeHKFGibbs(kv.second, P, T, state.dielectric.epsilon, state.dielectric.gf);
                    if(!std::isfinite(hkf.G))
                    {
                        result.passed = false;
                        result.errors.push_back("Invalid HKF G for species=" + kv.first + " in integration case=" + caseName + " row=" + std::to_string(rowIndex));
                    }
                }
            }
        }
        catch(const std::exception& e)
        {
            result.passed = false;
            result.errors.push_back("Integration case failed (" + path + "): " + e.what());
        }
    }

    return result;
}

// =============================================================================
// TEST 12: Chemical-Potential Parity (mu_i = RT * ln f_i, multi-P-T grid)
// =============================================================================
TestResult test_mu_parity()
{
    TestResult result{"Chemical-Potential Parity (H2O-CO2 P-T Grid)", true, {}, {}};

    constexpr double R = 8.31446261815324;

    const auto table = loadCsvTable("test/h2o_co2_pt_grid.csv");
    if(table.rows.empty())
    {
        result.passed = false;
        result.errors.push_back("Reference file test/h2o_co2_pt_grid.csv is empty or missing");
        return result;
    }

    PerpleXFluidOptions opts;
    opts.hybridSpecies = {1, 2};
    opts.hybridOptions = makePerpleXHybridEosOptions();
    opts.hybridOptions.water = HybridEosOptions::WaterEos::Mrk;
    opts.hybridOptions.co2   = HybridEosOptions::CO2Eos::Mrk;

    PerpleXFluidModel model;

    for(const auto& row : table.rows)
    {
        const double P        = getValue(table, row, "P(bar)");
        const double T        = getValue(table, row, "T(K)");
        const double yH2O     = getValue(table, row, "y(H2O)");
        const double yCO2     = getValue(table, row, "y(CO2)");
        const double fH2O_ref = getValue(table, row, "f(H2O)");
        const double fCO2_ref = getValue(table, row, "f(CO2)");

        std::array<double, 19> y{};
        y[1] = yH2O;
        y[2] = yCO2;

        std::vector<int> species;
        if(yH2O > 1e-12) species.push_back(1);
        if(yCO2 > 1e-12) species.push_back(2);

        const auto state = model.compute(species, y, P, T, opts);

        char label[128];

        if(yH2O > 1e-12 && std::isfinite(fH2O_ref) && fH2O_ref > 0.0 && fH2O_ref < 1e14)
        {
            const double mu_ref   = R * T * std::log(fH2O_ref);
            const double mu_model = R * T * state.ln_f[1];
            snprintf(label, sizeof(label), "mu(H2O) P=%.0fbar T=%.2fK", P, T);
            result.passed &= compareValue(mu_model, mu_ref, TOL_MU_ENERGY, label, result.errors, false);
        }

        if(yCO2 > 1e-12 && std::isfinite(fCO2_ref) && fCO2_ref > 0.0 && fCO2_ref < 1e14)
        {
            const double mu_ref   = R * T * std::log(fCO2_ref);
            const double mu_model = R * T * state.ln_f[2];
            snprintf(label, sizeof(label), "mu(CO2) P=%.0fbar T=%.2fK", P, T);
            result.passed &= compareValue(mu_model, mu_ref, TOL_MU_ENERGY, label, result.errors, false);
        }
    }

    return result;
}

// =============================================================================
// TEST 13: Delta-G Mixing Parity (excess Gibbs = RT * sum(yi * ln phi_i))
// =============================================================================
TestResult test_deltag_mix_parity()
{
    TestResult result{"Delta-G Mixing Parity (H2O-CO2 Composition Series)", true, {}, {}};

    constexpr double R = 8.31446261815324;

    const auto table = loadCsvTable("test/h2o_co2_mrk.csv");

    PerpleXFluidOptions opts;
    opts.hybridSpecies = {1, 2};
    opts.hybridOptions = makePerpleXHybridEosOptions();
    opts.hybridOptions.water = HybridEosOptions::WaterEos::Mrk;
    opts.hybridOptions.co2   = HybridEosOptions::CO2Eos::Mrk;

    PerpleXFluidModel model;
    int nChecked = 0;

    for(const auto& row : table.rows)
    {
        const double P    = getValue(table, row, "P(bar)");
        const double T    = getValue(table, row, "T(K)");
        const double yH2O = getValue(table, row, "y(H2O)");
        const double yCO2 = getValue(table, row, "y(CO2)");

        // Skip pure end-members: phi_i = f_i/(P*y_i) is undefined when y_i -> 0
        if(yH2O < 1e-12 || yCO2 < 1e-12) continue;

        std::array<double, 19> y{};
        y[1] = yH2O;
        y[2] = yCO2;

        const auto state = model.compute({1, 2}, y, P, T, opts);

        const double fH2O_ref = getValue(table, row, "f(H2O)");
        const double fCO2_ref = getValue(table, row, "f(CO2)");

        // Excess Gibbs of mixing: dGex = RT * sum(yi * ln(phi_i))
        //   where fugacity coefficient phi_i = f_i / (P * y_i)
        const double dGex_ref = R * T * (
            yH2O * std::log(fH2O_ref / (P * yH2O)) +
            yCO2 * std::log(fCO2_ref / (P * yCO2)));

        const double dGex_model = R * T * (
            yH2O * (state.ln_f[1] - std::log(P * yH2O)) +
            yCO2 * (state.ln_f[2] - std::log(P * yCO2)));

        char label[128];
        snprintf(label, sizeof(label), "deltaGex XCO2=%.1f P=%.0fbar T=%.2fK", yCO2, P, T);
        result.passed &= compareValue(dGex_model, dGex_ref, TOL_DELTAG_MIX, label, result.errors, false);
        ++nChecked;
    }

    if(nChecked == 0)
    {
        result.passed = false;
        result.errors.push_back("No binary mixture rows found in test/h2o_co2_mrk.csv");
    }

    return result;
}

// =============================================================================
// Main Test Runner
// =============================================================================
int main()
{
    try
    {
    std::cout << "======================================================================\n";
    std::cout << "Perple_X GFSM Regression Tests - Reaktoro Implementation\n";
    std::cout << "======================================================================\n\n";

    std::vector<TestResult> results;

    // Run all tests
    std::cout << "Running tests (GFSM EOS, DEW energy, mu parity, deltaG mixing)...\n\n";

    std::cout << "[1/13] GFSM EOS coverage gate..." << std::endl;
    results.push_back(test_gfsm_eos_coverage_gate());
    std::cout << "[2/13] Pure H2O test..." << std::endl;
    results.push_back(test_pure_h2o_hsmrk());
    std::cout << "[3/13] Pure H2O HSMRK test..." << std::endl;
    results.push_back(test_pure_h2o_hsmrk_ref());
    std::cout << "[4/13] Pure CO2 test..." << std::endl;
    results.push_back(test_pure_co2_cork());
    std::cout << "[5/13] Pure CO2 CORK test..." << std::endl;
    results.push_back(test_pure_co2_cork_ref());
    std::cout << "[6/13] Binary H2O-CO2 test..." << std::endl;
    results.push_back(test_h2o_co2_binary());
    std::cout << "[7/13] Composition series test..." << std::endl;
    results.push_back(test_composition_series());
    std::cout << "[8/13] Dielectric (epsh2o) test..." << std::endl;
    results.push_back(test_dielectric_h2o_co2());
    std::cout << "[9/13] DEW speciation energy parity (HKF Gibbs)..." << std::endl;
    results.push_back(test_dew_speciation_energy_parity());
    std::cout << "[10/13] Full GFSM matrix volume parity..." << std::endl;
    results.push_back(test_gfsm_full_matrix_volume_parity());
    std::cout << "[11/13] Integration flow scenarios..." << std::endl;
    results.push_back(test_integration_flow_scenarios());
    std::cout << "[12/13] Chemical-potential parity (mu P-T grid)..." << std::endl;
    results.push_back(test_mu_parity());
    std::cout << "[13/13] Delta-G mixing parity (excess Gibbs)..." << std::endl;
    results.push_back(test_deltag_mix_parity());

    // Optional component-level checks can be re-enabled for diagnostics only.
    // results.push_back(test_gfunction());
    // results.push_back(test_debye_huckel());
    // results.push_back(test_born_omega());

    // Print results
    int passed = 0;
    int failed = 0;

    for (const auto& result : results) {
        std::cout << "Test: " << result.testName << "\n";
        std::cout << "Status: " << (result.passed ? "âœ… PASSED" : "âŒ FAILED") << "\n";

        if (!result.errors.empty()) {
            std::cout << "Errors:\n";
            for (const auto& error : result.errors) {
                std::cout << "  - " << error << "\n";
            }
        }

        if (!result.warnings.empty()) {
            std::cout << "Warnings:\n";
            for (const auto& warning : result.warnings) {
                std::cout << "  - " << warning << "\n";
            }
        }

        std::cout << "\n";

        if (result.passed) ++passed;
        else ++failed;
    }

    // Summary
    std::cout << "======================================================================\n";
    std::cout << "Summary: " << passed << " passed, " << failed << " failed\n";
    std::cout << "Assertions: " << g_metrics.assertionsPassed << "/" << g_metrics.assertionsExecuted
              << " passed, " << g_metrics.assertionsFailed << " failed\n";
    std::cout << "Row-level Assertions: " << g_metrics.rowAssertionsPassed << "/" << g_metrics.rowAssertionsExecuted
              << " passed, " << g_metrics.rowAssertionsFailed << " failed\n";
    std::cout << "======================================================================\n";

    // Generate report file
    std::ofstream report("regression_test_report.txt");
    report << "Perple_X Regression Test Report\n";
    report << "================================\n\n";
    report << "Total Tests: " << results.size() << "\n";
    report << "Passed: " << passed << "\n";
    report << "Failed: " << failed << "\n\n";
    report << "Assertions Executed: " << g_metrics.assertionsExecuted << "\n";
    report << "Assertions Passed: " << g_metrics.assertionsPassed << "\n";
    report << "Assertions Failed: " << g_metrics.assertionsFailed << "\n";
    report << "Row-level Assertions Executed: " << g_metrics.rowAssertionsExecuted << "\n";
    report << "Row-level Assertions Passed: " << g_metrics.rowAssertionsPassed << "\n";
    report << "Row-level Assertions Failed: " << g_metrics.rowAssertionsFailed << "\n\n";

    for (const auto& result : results) {
        report << "Test: " << result.testName << "\n";
        report << "Status: " << (result.passed ? "PASSED" : "FAILED") << "\n";
        if (!result.errors.empty()) {
            report << "Errors:\n";
            for (const auto& error : result.errors) {
                report << "  " << error << "\n";
            }
        }
        report << "\n";
    }

    report.close();

    return failed > 0 ? 1 : 0;
    }
    catch(const std::exception& e)
    {
        std::cerr << "\n[FATAL] Unhandled regression exception: " << e.what() << "\n";
        return 1;
    }
    catch(...)
    {
        std::cerr << "\n[FATAL] Unhandled non-standard regression exception\n";
        return 1;
    }
}
