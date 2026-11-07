#Requires -Version 7

<#
.SYNOPSIS
  Build (in WSL) and run the Brane dev container on a Windows host with Podman.

.DESCRIPTION
  Assumption: `mise`, `wsl`, and `podman` are already installed and working on
  the host — this script does not install anything.

  The pinned Lean 4 toolchain is Linux-only, so the OCI image is built by
  `mise oci build` inside WSL (see the root mise.toml [oci] section). This
  script drives the whole flow:

    1. `mise oci build` in WSL, writing the OCI layout to <RepoDir>/mise-oci
       (skipped if the layout already exists, unless -Rebuild).
    2. Pack the layout as an OCI archive tar and `podman load` it.
    3. `podman run` an interactive bash in the image with the repo mounted at
       /workspace (or run -Command non-interactively).

  Once inside the container, for example:
     cd /workspace/projects/metropolis && lake build && lake exe metropolis

.PARAMETER RepoDir
  The brane repository checkout. Defaults to the parent of this script's
  directory, which is the repo root in a normal clone.

.PARAMETER Rebuild
  Force a fresh `mise oci build` even if mise-oci/ already has an image.

.PARAMETER WslDistro
  Specific WSL distribution to use (default: the WSL default distro).

.PARAMETER Command
  If given, run this command in the container instead of an interactive bash.

.EXAMPLE
  .\scripts\start-brane.ps1
  .\scripts\start-brane.ps1 -Rebuild
  .\scripts\start-brane.ps1 -Command "lake build && lake exe metropolis"
#>

param(
  [string]$RepoDir = (Split-Path -Parent $PSScriptRoot),
  [switch]$Rebuild,
  [string]$WslDistro = "",
  [string]$Command = ""
)

$ErrorActionPreference = "Stop"

# --- Resolve & validate paths -------------------------------------------------
$RepoRoot = $RepoDir
$ImageDirWin = Join-Path $RepoRoot "mise-oci"
$TarPathWin = Join-Path $env:TEMP ("brane-oci-{0}.tar" -f (Get-Date -Format yyyyMMddHHmmss))

if (-not (Test-Path (Join-Path $RepoRoot "mise.toml"))) {
  throw "Bad repo dir: $RepoRoot (no mise.toml found)"
}

## utility: run a shell line inside WSL, with Windows path translated first
function Invoke-Wsl {
  param([Parameter(Mandatory)][string]$WinPath, [Parameter(Mandatory)][string]$ShellCmd)
  $wslArgs = @()
  if ($WslDistro) { $wslArgs = @("--distribution", $WslDistro) }
  $wslArgs += @("--exec", "bash", "-lc", ('cd "$(wslpath -a ' + "'" + $WinPath + "'" + ')" && ' + $ShellCmd))
  & wsl @wslArgs
  if ($LASTEXITCODE -ne 0) { throw "WSL command failed (exit $LASTEXITCODE): $ShellCmd" }
}

Write-Host "== Repo      : $RepoRoot"
Write-Host "== Image dir : $ImageDirWin"

# --- 1. Build the image in WSL (Linux toolchain binaries required) -------------
$needBuild = $Rebuild -or (-not (Test-Path (Join-Path $ImageDirWin "index.json")))
if ($needBuild) {
  Write-Host "== Building image via mise oci (WSL) ..."
  Invoke-Wsl -WinPath $RepoRoot -ShellCmd "mise install && mise oci build -o ./mise-oci"
} else {
  Write-Host "== mise-oci/ layout already exists — skipping build (use -Rebuild to force)"
}

# --- 2. Pack the OCI layout as a tar and load it into Podman -------------------
Write-Host "== Packing layout into OCI archive: $TarPathWin"
Invoke-Wsl -WinPath $RepoRoot `
  -ShellCmd ('tar -cf "$(wslpath -a ' + "'" + $TarPathWin + "'" + ')" -C ./mise-oci .')

Write-Host "== Loading image into podman ..."
podman load -i $TarPathWin
if ($LASTEXITCODE -ne 0) { throw "podman load failed (exit $LASTEXITCODE)" }
Remove-Item $TarPathWin -ErrorAction SilentlyContinue

# --- 3. Run --------------------------------------------------------------------
# The OCI archive carries the ref name from mise.toml [oci].tag.
$imageRef = "ghcr.io/TreeHappy/brane:latest"
$repoMount = "${RepoRoot}:/workspace"

if ($Command) {
  Write-Host "== Running in container: $Command"
  & podman run --rm `
      --volume $repoMount `
      --workdir /workspace `
      $imageRef `
      bash -lc $Command
} else {
  Write-Host "== Interactive shell in $imageRef (repo at /workspace)"
  & podman run --rm -it `
      --volume $repoMount `
      --workdir /workspace `
      $imageRef `
      bash
}
if ($LASTEXITCODE -ne 0) { throw "podman run failed (exit $LASTEXITCODE)" }
