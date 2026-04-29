[CmdletBinding()]
param(
    [ValidateSet('standard', 'dew')]
    [string]$Mode = 'standard',

    [string]$BuildDir,

    [ValidateSet('Debug', 'Release', 'RelWithDebInfo', 'MinSizeRel')]
    [string]$Config = 'Release',

    [string]$Generator = 'Visual Studio 16 2019',

    [string]$Architecture = 'x64',

    [string]$InstallPrefix,

    [switch]$SkipEnvironment,

    [switch]$SkipInstall,

    [switch]$SkipTests,

    [switch]$RunRegressionSuite,

    [string[]]$RegressionTags = @('smoke')
)

$ErrorActionPreference = 'Stop'
$repoRoot = Split-Path -Parent $PSScriptRoot

if (-not $BuildDir) {
    $BuildDir = 'build'
}

function Get-CondaExecutable {
    if ($env:CONDA_EXE -and (Test-Path $env:CONDA_EXE)) {
        return $env:CONDA_EXE
    }

    try {
        return (Get-Command conda -ErrorAction Stop).Source
    }
    catch {
        throw 'Could not find conda. Install Miniconda/Anaconda and ensure conda is available in PATH or CONDA_EXE is set.'
    }
}

function Invoke-External {
    param(
        [Parameter(Mandatory = $true)]
        [string]$FilePath,

        [Parameter(Mandatory = $true)]
        [string[]]$Arguments,

        [Parameter(Mandatory = $true)]
        [string]$Description
    )

    Write-Host "==> $Description" -ForegroundColor Cyan
    & $FilePath @Arguments

    if ($LASTEXITCODE -ne 0) {
        throw "$Description failed with exit code $LASTEXITCODE."
    }
}

function Invoke-Conda {
    param(
        [Parameter(Mandatory = $true)]
        [string[]]$Arguments,

        [Parameter(Mandatory = $true)]
        [string]$Description
    )

    Invoke-External -FilePath $script:CondaExe -Arguments $Arguments -Description $Description
}

function Invoke-ReaktoroEnv {
    param(
        [Parameter(Mandatory = $true)]
        [string[]]$Arguments,

        [Parameter(Mandatory = $true)]
        [string]$Description
    )

    $fullArguments = @('run', '--no-capture-output', '-n', 'reaktoro') + $Arguments
    Invoke-Conda -Arguments $fullArguments -Description $Description
}

$script:CondaExe = Get-CondaExecutable

$vsWhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
if (-not (Test-Path $vsWhere)) {
    Write-Warning 'vswhere.exe was not found. Install Visual Studio 2019 Build Tools with the C++ workload if CMake cannot find the generator.'
}

Push-Location $repoRoot
try {
    if (-not $SkipEnvironment) {
        Invoke-Conda -Arguments @('install', '-n', 'base', '-c', 'conda-forge', '-y', 'mamba', 'conda-devenv') -Description 'Installing mamba and conda-devenv in the base environment'
        Invoke-Conda -Arguments @('devenv', '-e', 'mamba', '-f', 'environment.devenv.yml') -Description 'Creating or updating the reaktoro Conda environment'
    }

    $cmakeConfigureArgs = @(
        'cmake',
        '-S', '.',
        '-B', $BuildDir,
        '-G', $Generator,
        '-A', $Architecture,
        '-DREAKTORO_BUILD_PYTHON=ON',
        '-DREAKTORO_BUILD_TESTS=ON',
        '-DREAKTORO_BUILD_DOCS=OFF'
    )

    if ($InstallPrefix) {
        $cmakeConfigureArgs += "-DCMAKE_INSTALL_PREFIX=$InstallPrefix"
    }

    Invoke-ReaktoroEnv -Arguments $cmakeConfigureArgs -Description "Configuring the $BuildDir build directory"
    Invoke-ReaktoroEnv -Arguments @('cmake', '--build', $BuildDir, '--config', $Config, '--parallel') -Description "Building Reaktoro in $BuildDir"

    if (-not $SkipInstall) {
        Invoke-ReaktoroEnv -Arguments @('cmake', '--install', $BuildDir, '--config', $Config) -Description "Installing Reaktoro from $BuildDir"
    }

    if (-not $SkipTests) {
        Invoke-ReaktoroEnv -Arguments @('ctest', '--test-dir', $BuildDir, '-C', $Config, '--output-on-failure') -Description "Running CTest from $BuildDir"
    }

    if ($RunRegressionSuite) {
        $regressionScript = Join-Path 'DEW_Experimental_Benchmark' 'regression_suite.py'
        $regressionArgs = @('python', $regressionScript, '--continue-on-fail')
        if ($RegressionTags -and $RegressionTags.Count -gt 0) {
            $regressionArgs += '--tags'
            $regressionArgs += $RegressionTags
        }
        Invoke-ReaktoroEnv -Arguments $regressionArgs -Description "Running DEW/PerplexDEW regression suite"
    }
}
finally {
    Pop-Location
}

Write-Host ''
Write-Host 'Setup completed.' -ForegroundColor Green
Write-Host "Build directory: $BuildDir"
Write-Host "Configuration: $Config"
Write-Host "Mode: $Mode"

if ($Mode -eq 'dew') {
    Write-Host 'The DEW benchmark scripts should now be able to find the local build in build.'
}
