param(
    [switch]$Clean,
    [switch]$SkipRun,
    [string]$Configuration = "Release",
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$SimulatorArgs
)

$ErrorActionPreference = "Stop"

$Root = $PSScriptRoot
$Sources = Get-ChildItem -Path (Join-Path $Root "src") -Filter "*.cpp" | Select-Object -ExpandProperty FullName
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

if ($Sources.Count -eq 0) {
    throw "No C++ source files found in src/"
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
        "-flto",
        "-std=c++17",
        "-fopenmp",
        "-Wall",
        "-Wextra",
        "-o", $Exe
    )
    $args += $Sources
    $args += "-lm"
    Invoke-Checked $gpp $args
}
elseif ($clangpp) {
    $compilerBin = Split-Path -Parent $clangpp
    if ($env:Path -notlike "*$compilerBin*") {
        $env:Path = "$compilerBin;$env:Path"
    }
    $args = @(
        "-O2",
        "-flto",
        "-std=c++17",
        "-Wall",
        "-Wextra",
        "-o", $Exe
    )
    $args += $Sources
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
        "/GL",
        "/openmp",
        "/Fe:$Exe",
        "/Fo:$objDir\"
    )
    $args += $Sources
    Invoke-Checked $cl $args
}
else {
    throw "No supported C++ compiler found on PATH. Install g++, clang++, or run from a Visual Studio Developer PowerShell with cl.exe available."
}

Write-Host "Built: $Exe"
Unblock-File -Path $Exe -ErrorAction SilentlyContinue
Start-Sleep -Seconds 1

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
