#!/usr/bin/env bash
# Build the in-process Thermo RAW helper (vulcan/io/dotnet/VulcanRaw) into tools/vulcanraw/.
# Installs the .NET 8 SDK into tools/dotnet-sdk/ if missing (dotnet-install.sh, no root needed).
# References RawFileReader DLLs shipped with ThermoRawFileParser in conda_envs/trfp (Thermo licence).
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SDK="$ROOT/tools/dotnet-sdk"
export DOTNET_SYSTEM_GLOBALIZATION_INVARIANT=1 DOTNET_CLI_TELEMETRY_OPTOUT=1 DOTNET_NOLOGO=1 DOTNET_ROOT="$SDK"
export DOTNET_CLI_HOME="${TMPDIR:-/tmp}/vulcan_dotnet_home" NUGET_PACKAGES="${TMPDIR:-/tmp}/vulcan_nuget"
if [ ! -x "$SDK/dotnet" ]; then
  curl -sSL -o "${TMPDIR:-/tmp}/dotnet-install.sh" https://dot.net/v1/dotnet-install.sh
  bash "${TMPDIR:-/tmp}/dotnet-install.sh" --channel 8.0 --install-dir "$SDK" --no-path
fi
cd "$ROOT/vulcan/io/dotnet/VulcanRaw"
"$SDK/dotnet" build -c Release -o "$ROOT/tools/vulcanraw" -p:TrfpBin="${VULCAN_TRFP_BIN:-$ROOT/conda_envs/trfp/bin}"
rm -rf obj bin
echo "built $ROOT/tools/vulcanraw/VulcanRaw.dll"
