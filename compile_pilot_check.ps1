# Direct compile check for ActivityModelMAGEMinSolidSolutionPilot.cpp
$DEW_PATH = 'C:/Users/stanroozen/Documents/Projects/reaktoro-dev/reaktoro/embedded/databases/DEW/dew2024-aqueous.yaml'
$SRC = 'C:/Users/stanroozen/Documents/Projects/reaktoro-dev/reaktoro/Reaktoro/Models/ActivityModels/ActivityModelMAGEMinSolidSolutionPilot.cpp'
$INCLUDE_BASE = 'C:/Users/stanroozen/Documents/Projects/reaktoro-dev/reaktoro'
$INCLUDE_BUILD = 'C:/Users/stanroozen/Documents/Projects/reaktoro-dev/reaktoro/build/_cmrc/include'
$SYSTEM_INCLUDE = 'C:/Users/stanroozen/anaconda3/envs/reaktoro/Library/include'
$EIGEN_INCLUDE = 'C:/Users/stanroozen/anaconda3/envs/reaktoro/Library/include/eigen3'

$args_list = @(
    "-DAUTODIFF_ENABLE_IMPLICIT_CONVERSION_REAL=1",
    "-DDEW_AQUEOUS_DB_PATH=`"$DEW_PATH`"",
    "-DFMT_SHARED",
    "-DReaktoro_EXPORTS",
    "-DSPDLOG_COMPILED_LIB",
    "-DSPDLOG_FMT_EXTERNAL",
    "-DSPDLOG_SHARED_LIB",
    "-DUSE_SPDLOG_PRECOMPILED",
    "-I$INCLUDE_BASE",
    "-I$INCLUDE_BUILD",
    "-isystem", $SYSTEM_INCLUDE,
    "-isystem", $EIGEN_INCLUDE,
    "-O0",
    "-DNDEBUG",
    "-std=c++20",
    "-fsyntax-only",
    $SRC
)

Write-Host "Running compiler..."
$proc = Start-Process -FilePath "C:\msys64\ucrt64\bin\c++.exe" `
    -ArgumentList $args_list `
    -NoNewWindow -Wait -PassThru `
    -RedirectStandardError "c:/Users/stanroozen/Documents/Projects/reaktoro-dev/reaktoro/compile_errors.txt"

Write-Host "Exit code: $($proc.ExitCode)"
Get-Content "c:/Users/stanroozen/Documents/Projects/reaktoro-dev/reaktoro/compile_errors.txt" | Select-Object -First 60
