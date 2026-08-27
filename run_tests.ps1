#!/usr/bin/env powershell
<#
.SYNOPSIS
    Unified test runner for Reaktoro test suite.

.DESCRIPTION
    Discovers Python environment (3.12), auto-detects reaktoro4py.pyd location,
    and runs pytest with appropriate markers and options.

.PARAMETER Tier
    Test execution tier:
    - 1 (default): Workflow coverage — 18 critical tests (~5s)
    - 2: Extended smoke tests — 30 tests (~60s)
    - 3: Full test stack — all tests including unit/bindings (~120s)
    - specific: Run a specific test file or test name

.PARAMETER Filter
    Filter tests by name (passed to pytest -k)
    Example: -Filter "dew" runs tests matching "dew"

.PARAMETER Verbose
    Show full pytest output (verbose mode)

.EXAMPLE
    # Run workflow coverage (default)
    .\run_tests.ps1

    # Run all smoke tests
    .\run_tests.ps1 -Tier 2

    # Run GFSM handoff tests
    .\run_tests.ps1 -Filter gfsm

    # Run specific test
    .\run_tests.ps1 -Filter "test_gfsm_handoff_consumed_by_perplexdew"

    # Verbose output
    .\run_tests.ps1 -Tier 2 -Verbose
#>

param(
    [ValidateSet('1', '2', '3', 'specific')]
    [string]$Tier = '1',

    [string]$Filter = '',
    [switch]$Verbose
)

# Find Python 3.12 environment
$python_exe = $null
$envs = @(
    'C:\Users\stanroozen\anaconda3\envs\reaktoro\python.exe',
    'C:\Users\stanroozen\anaconda3\python.exe',
    $(python -c "import sys; print(sys.executable)" 2>$null)
)

foreach ($py in $envs) {
    if (Test-Path $py) {
        $version = & $py --version 2>&1
        if ($version -match 'Python 3\.12') {
            $python_exe = $py
            break
        }
    }
}

if (-not $python_exe) {
    Write-Host "ERROR: Python 3.12 not found!" -ForegroundColor Red
    Write-Host "Please activate the 'reaktoro' conda environment:" -ForegroundColor Yellow
    Write-Host "  conda activate reaktoro" -ForegroundColor Cyan
    exit 1
}

Write-Host "Using Python: $python_exe" -ForegroundColor Green
& $python_exe --version

# Set up test command based on tier
$repo_root = Split-Path -Parent $PSScriptRoot
$test_cmd = @($python_exe, '-m', 'pytest')

if ($Verbose) {
    $test_cmd += '-vv'
} else {
    $test_cmd += '-v'
}

switch ($Tier) {
    '1' {
        Write-Host "`nRunning Tier 1: Workflow Coverage (18 tests, ~5s)" -ForegroundColor Cyan
        $test_cmd += @('Testing/', '-m', 'workflow_coverage')
    }
    '2' {
        Write-Host "`nRunning Tier 2: Extended Smoke Tests (30 tests, ~60s)" -ForegroundColor Cyan
        $test_cmd += @('Testing/regression/smoke/', '--ignore=Testing/scripts', '--ignore=Testing/unit', '--ignore=Testing/bindings')
    }
    '3' {
        Write-Host "`nRunning Tier 3: Full Test Stack (~120s)" -ForegroundColor Cyan
        $test_cmd += @('Testing/', '--ignore=Testing/scripts')
    }
    'specific' {
        if (-not $Filter) {
            Write-Host "ERROR: -Filter required for -Tier specific" -ForegroundColor Red
            exit 1
        }
        Write-Host "`nRunning specific tests matching: $Filter" -ForegroundColor Cyan
        $test_cmd += @('Testing/', '-k', $Filter)
    }
}

# Apply filter if provided
if ($Filter -and $Tier -ne 'specific') {
    $test_cmd += @('-k', $Filter)
}

# Add options
$test_cmd += @('--tb=short', '--timeout=300')

# Run tests
Write-Host "`nCommand: $($test_cmd -join ' ')" -ForegroundColor Yellow
Write-Host ""

& $python_exe $($test_cmd | Select-Object -Skip 1)
$exit_code = $LASTEXITCODE

# Summary
Write-Host ""
if ($exit_code -eq 0) {
    Write-Host "✓ Tests passed!" -ForegroundColor Green
} else {
    Write-Host "✗ Tests failed (exit code: $exit_code)" -ForegroundColor Red
}

exit $exit_code
