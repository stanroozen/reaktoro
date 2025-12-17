#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <cmath>
#include <limits>

// ------------------------- Small CSV helper -------------------------

struct CsvRow {
    std::vector<std::string> fields;
};

inline void trim_inplace(std::string& s)
{
    auto is_space = [](unsigned char c) { return std::isspace(c) != 0; };

    std::size_t start = 0;
    while (start < s.size() && is_space(static_cast<unsigned char>(s[start])))
        ++start;

    std::size_t end = s.size();
    while (end > start && is_space(static_cast<unsigned char>(s[end - 1])))
        --end;

    if (start == 0 && end == s.size())
        return;

    s = s.substr(start, end - start);
}

// CSV loader with proper quoted field handling for fields containing commas
inline std::vector<CsvRow> load_csv(const std::string& path, bool skip_header)
{
    std::ifstream in(path);
    if (!in)
        throw std::runtime_error("Cannot open CSV file: " + path);

    std::vector<CsvRow> rows;
    std::string line;

    bool is_first = true;
    while (std::getline(in, line)) {
        if (line.empty())
            continue;
        if (is_first && skip_header) {
            is_first = false;
            continue;
        }
        is_first = false;

        CsvRow row;
        std::string cell;
        bool in_quotes = false;

        for (size_t i = 0; i < line.size(); ++i) {
            char c = line[i];

            if (c == '"') {
                in_quotes = !in_quotes;
                cell += c; // Keep quotes for strip_quotes to remove later
            }
            else if (c == ',' && !in_quotes) {
                trim_inplace(cell);
                row.fields.push_back(cell);
                cell.clear();
            }
            else {
                cell += c;
            }
        }

        // Don't forget the last field
        trim_inplace(cell);
        row.fields.push_back(cell);

        rows.push_back(std::move(row));
    }
    return rows;
}

// Remove surrounding quotes from a string field (e.g., "ACETATE,AQ" -> ACETATE,AQ)
inline std::string strip_quotes(const std::string& s)
{
    std::string result = s;
    trim_inplace(result);

    // Remove surrounding double quotes
    if (result.size() >= 2 && result.front() == '"' && result.back() == '"') {
        result = result.substr(1, result.size() - 2);
    }

    return result;
}

// Parse a double if the field is non-empty and not NaN.
// Returns false if the field should be treated as "missing".
inline bool parse_maybe_double(const std::string& s, double& value)
{
    std::string t = s;
    trim_inplace(t);
    if (t.empty())
        return false;

    // Handle NaN-like markers
    for (char& c : t) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (t == "nan" || t == "NaN" || t == "nan(ind)" || t == "inf" || t == "-inf")
        return false;

    try {
        value = std::stod(t);
    } catch (...) {
        return false;
    }
    return true;
}

// Simple rel+abs tolerance comparison
inline bool almost_equal(double a, double b, double abs_tol, double rel_tol)
{
    double diff = std::fabs(a - b);
    if (diff <= abs_tol)
        return true;

    double scale = std::max(std::fabs(a), std::fabs(b));
    if (scale == 0.0)
        return diff <= abs_tol;

    if (diff <= rel_tol * scale)
        return true;

    return false;
}

