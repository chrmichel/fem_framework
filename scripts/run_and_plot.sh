#!/usr/bin/env bash
set -euo pipefail

# Adjust if your build dir or target name differs
BUILD_DIR="${BUILD_DIR:-build}"
APP_PATH="${APP_PATH:-$BUILD_DIR/apps/poisson1d/poisson1d}"

if [ ! -x "$APP_PATH" ]; then
  echo "App not found or not executable: $APP_PATH"
  echo "Build first: cmake --build $BUILD_DIR"
  echo "Or set APP_PATH explicitly, e.g.: APP_PATH=build/apps/... ./scripts/run_and_plot.sh"
  exit 1
fi

# Run app from repo root so results/ lands in repo root (recommended)
"$APP_PATH"

# Ensure results exists
if [ ! -d "results" ]; then
  echo "No results/ directory found after running app."
  exit 1
fi

# Plot latest results
gnuplot -persist scripts/plot_latest.gp