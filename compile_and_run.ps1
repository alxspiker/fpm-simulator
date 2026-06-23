param(
    [switch]$Clean,
    [switch]$SkipRun,
    [string]$Configuration = "Release",
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$SimulatorArgs
)

$ErrorActionPreference = "Stop"

$Root = $PSScriptRoot
$Source = Join-Path $Root "src\fpm_axcore_simulator.cpp"
$BuildDir = Join-Path $Root "build"
$Exe = Join-Path $BuildDir "fpm_axcore.exe"

function Invoke-Checked {
    param(
        [string]$FilePath,
        [string[]]$Arguments
    )

    Write-Host "> $FilePath $($Arguments -join ' ')"
    & $FilePath @Arguments
    if ($LASTEXITCODE -ne 0) {
        throw "Command failed with exit code $LASTEXITCODE."
    }
}

function Find-Tool {
    param(
        [string]$Name,
        [string[]]$FallbackPaths
    )

    $cmd = Get-Command $Name -ErrorAction SilentlyContinue
    if ($cmd) {
        return $cmd.Source
    }

    foreach ($path in $FallbackPaths) {
        if (Test-Path $path) {
            return $path
        }
    }

    return $null
}

if (!(Test-Path $Source)) {
    throw "Source file not found: $Source"
}

if ($Clean -and (Test-Path $BuildDir)) {
    Remove-Item -LiteralPath $BuildDir -Recurse -Force
}

New-Item -ItemType Directory -Force -Path $BuildDir | Out-Null

$gpp = Find-Tool "g++" @(
    "C:\msys64\ucrt64\bin\g++.exe",
    "C:\msys64\mingw64\bin\g++.exe",
    "C:\msys64\clang64\bin\g++.exe",
    "C:\mingw64\bin\g++.exe"
)
$clangpp = Find-Tool "clang++" @(
    "C:\msys64\clang64\bin\clang++.exe",
    "C:\msys64\ucrt64\bin\clang++.exe",
    "C:\Program Files\LLVM\bin\clang++.exe"
)
$cl = Find-Tool "cl.exe" @()

Write-Host "Building FPM AxCore simulator ($Configuration)..."

if ($gpp) {
    $compilerBin = Split-Path -Parent $gpp
    if ($env:Path -notlike "*$compilerBin*") {
        $env:Path = "$compilerBin;$env:Path"
    }
    $args = @(
        "-O2",
        "-std=c++17",
        "-fopenmp",
        "-Wall",
        "-Wextra",
        "-o", $Exe,
        $Source,
        "-lm"
    )
    Invoke-Checked $gpp $args
}
elseif ($clangpp) {
    $compilerBin = Split-Path -Parent $clangpp
    if ($env:Path -notlike "*$compilerBin*") {
        $env:Path = "$compilerBin;$env:Path"
    }
    $args = @(
        "-O2",
        "-std=c++17",
        "-Wall",
        "-Wextra",
        "-o", $Exe,
        $Source
    )
    Invoke-Checked $clangpp $args
}
elseif ($cl) {
    $objDir = Join-Path $BuildDir "obj"
    New-Item -ItemType Directory -Force -Path $objDir | Out-Null
    $args = @(
        "/nologo",
        "/EHsc",
        "/std:c++17",
        "/O2",
        "/openmp",
        "/Fe:$Exe",
        "/Fo:$objDir\",
        $Source
    )
    Invoke-Checked $cl $args
}
else {
    throw "No supported C++ compiler found on PATH. Install g++, clang++, or run from a Visual Studio Developer PowerShell with cl.exe available."
}

Write-Host "Built: $Exe"

if (!$SkipRun) {
    Write-Host "Running simulator..."
    Push-Location $Root
    try {
        Invoke-Checked $Exe $SimulatorArgs
    }
    finally {
        Pop-Location
    }
}
