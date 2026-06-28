param(
    [switch]$Clean,
    [switch]$RunSmoke,
    [string]$Optimization = "-O0",
    [string]$EmsdkPath = "$env:USERPROFILE\emsdk",
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$SimulatorArgs
)

$ErrorActionPreference = "Stop"

$Root = $PSScriptRoot
$Sources = Get-ChildItem -Path (Join-Path $Root "src") -Filter "*.cpp" | Select-Object -ExpandProperty FullName
$OutDir = Join-Path $Root "web"
$OutJs = Join-Path $OutDir "fpm_axcore.js"
$OutWasm = Join-Path $OutDir "fpm_axcore.wasm"
$EmsdkEnv = Join-Path $EmsdkPath "emsdk_env.bat"

if ($Sources.Count -eq 0) {
    throw "No C++ source files found in src/"
}

if (!(Test-Path $EmsdkEnv)) {
    throw "Emscripten SDK not found at $EmsdkPath. Install with: git clone https://github.com/emscripten-core/emsdk.git $EmsdkPath; cd $EmsdkPath; .\emsdk.bat install latest; .\emsdk.bat activate latest"
}

if ($Clean -and (Test-Path $OutDir)) {
    Remove-Item -LiteralPath $OutDir -Recurse -Force
}

New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

Write-Host "Building FPM AxCore simulator for WebAssembly..."

$compileArgs = @(
    $Optimization,
    "-flto",
    "-std=c++17",
    "-Wall",
    "-Wextra",
    "-s", "WASM=1",
    "-s", "ALLOW_MEMORY_GROWTH=1",
    "-s", "STACK_SIZE=8388608",
    "-s", "EXIT_RUNTIME=1",
    "-s", "NODERAWFS=1",
    "-s", "ENVIRONMENT=web,node",
    "-o", $OutJs
)
$compileArgs += $Sources

$quotedArgs = ($compileArgs | ForEach-Object { '"' + ($_ -replace '"', '\"') + '"' }) -join " "
$cmd = "`"$EmsdkEnv`" >nul && em++ $quotedArgs"
Write-Host "> em++ $($compileArgs -join ' ')"
cmd /c $cmd
if ($LASTEXITCODE -ne 0) {
    throw "WASM compile failed with exit code $LASTEXITCODE."
}

if (!(Test-Path $OutJs) -or !(Test-Path $OutWasm)) {
    throw "Expected WASM outputs were not produced: $OutJs / $OutWasm"
}

Write-Host "Built:"
Write-Host "  $OutJs"
Write-Host "  $OutWasm"

if ($RunSmoke) {
    $SmokeArgs = $SimulatorArgs
    if (!$SmokeArgs -or $SmokeArgs.Count -eq 0) {
        $SmokeArgs = @("--torsion-phase-lock-output", "web\torsion_smoke.json")
    }

    $nodeArgs = @($OutJs) + $SmokeArgs
    $quotedNodeArgs = ($nodeArgs | ForEach-Object { '"' + ($_ -replace '"', '\"') + '"' }) -join " "
    Write-Host "Running WASM smoke test..."
    Write-Host "> node $($nodeArgs -join ' ')"
    cmd /c "`"$EmsdkEnv`" >nul && node $quotedNodeArgs"
    if ($LASTEXITCODE -ne 0) {
        throw "WASM smoke test failed with exit code $LASTEXITCODE."
    }
}
